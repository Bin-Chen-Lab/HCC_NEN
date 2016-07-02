#compare patient tumor samples vs their drug response data
#whole-genome data will be used
#look proportion of patients who will be reversed by NEN/Niclo in vivo and in vitro.
#only look patients who have matched non-tumor

setwd("/Users/binchen1/Documents/stanford/hcc/data/")

library(DESeq) 
library(edgeR) 

cancer = "LIHC"
load(paste(cancer, "/tumor_info.RData", sep=""))
load(paste(cancer, "/countTable.RData", sep="")) #use scale esimate below
load(paste(cancer, "/cds.RData", sep=""))

#load(paste(cancer, "/scale_estimate.RData", sep=""))

data2 = read.csv("microarray//Group3/data2.csv")
data2 = subset(data2, fail.count<3)
drug_fc = aggregate(Test.Control.fc ~ GENE_ID, data2, mean)
names(drug_fc) = c( "GeneID", "fc")

#normalized_cds_all = countTable
normalized_cds_all = counts( cds, normalized=TRUE )

tumor_info = subset(tumor_info, quality != "outlier" & rownames(tumor_info) %in% colnames(normalized_cds_all))
gene_ids = sapply(rownames(normalized_cds_all), function(id){unlist(strsplit(id, "\\|"))[2]})


signatures <- read.table("LIHC/drug/LIHC_dz_signature_cmap_final_v1.txt",header=T,sep="\t") #at least cover gene id and value that can be either fold change or p value
dz_sig = subset(signatures, select=c("GeneID", "value"))

patient_freq = table(tumor_info$bcr_patient_barcode)
matched_patients = names(patient_freq[patient_freq>1])
#non tumor gene expression

patient_sigs = data.frame(gene_ids)
patient_drug_cor = data.frame()
for (patient in matched_patients){
   patient_tumor = rownames(tumor_info[tumor_info$bcr_patient_barcode == patient & tumor_info$type == "tumor",])
   tumor_sig = normalized_cds_all[, patient_tumor[1]]
   
   patient_non_tumor = rownames(tumor_info[tumor_info$bcr_patient_barcode == patient & tumor_info$type == "non-tumor",])
   if (length(patient_non_tumor)<1) { next }
   
   non_tumor_sig = normalized_cds_all[, patient_non_tumor]
   
   #compare to non tumors
   patient_sig = log(tumor_sig+1) - log(non_tumor_sig + 1)
   patient_sigs = cbind(patient_sigs, data.frame(patient_sig))
   
   patient_sig = data.frame(patient_sig, gene_ids)

   
   patient_drug = merge( drug_fc, patient_sig, by.x="GeneID", by.y="gene_ids")
   patient_drug = merge( dz_sig,patient_drug, by="GeneID", sort=F)
   
   cor = cor.test(patient_drug$fc, patient_drug$patient_sig, method="pearson")
   patient_drug_cor = rbind(patient_drug_cor, data.frame(patient, p=cor$p.value, cor=cor$estimate))
}
colnames(patient_sigs) = c("GeneID", matched_patients)

patient_drug_cor$p.adj = p.adjust(patient_drug_cor$p, "fdr")

sum(patient_drug_cor$p.adj<0.01)/nrow(patient_drug_cor)
#81%
save(patient_sigs, file="tcga_matched_patient_expr.RData")

patient_drug_cor = patient_drug_cor[order(patient_drug_cor$p.adj), ]
#clinical features of non-reversed patient
non_responder = tumor_info[tumor_info$bcr_patient_barcode %in%  as.character(patient_drug_cor$patient[patient_drug_cor$p.adj>0.01]),]
#View(non_responder)

patients_drug = merge( drug_fc, patient_sigs, by.x="GeneID", by.y="GeneID")
patients_drug = merge( dz_sig,patients_drug, by="GeneID", sort=F)


patients_drug_rank = patients_drug[,-c(1,2)]
for (i in 1:ncol(patients_drug_rank)){
  patients_drug_rank[, i ] = rank(-1 * patients_drug_rank[, i ] ) #upregulated ranked better
}

col_sorted = sort(cor(patients_drug[,-c(1,2)], method="pearson")["fc",-1])    
patients_drug_rank = patients_drug_rank[,c("fc", names(col_sorted))]
patients_drug_rank = patients_drug_rank[order(patients_drug_rank$fc), ]

#label patients by significance
patient_labels = sapply(names(col_sorted), function(patient){
  if (patient %in% non_responder$bcr_patient_barcode){
    patient
  }else{
    paste("*", patient, sep="")
  }
})

library(gplots)
library(ggplot2)

pdf("LIHC/patient_response.pdf")
colPal <- greenred(100)
par(mar=c(12, 4, 2, 0.5))
#image(t(druggable_targets_pos), col=redblue(2))
image(t(patients_drug_rank), col= colPal,   axes=F, srt=45)
axis(1,  at=seq(0,1,length.out= ncol( patients_drug_rank ) ), labels= F)
text(x = seq(0,1,length.out=ncol( patients_drug_rank ) ), c(-0.05),
     labels = c( "NEN",patient_labels), srt = 45, pos=2,offset=0.05, xpd = TRUE, cex=0.7)

dev.off()
