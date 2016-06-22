#find patients who may have efficacy to niclosamide

setwd("~/Documents/stanford/hcc/data/")

library(DESeq) 
library(edgeR) 
library(survival)
library(KMsurv)

cancer = "LIHC"
load(paste(cancer, "/tumor_info.RData", sep=""))

load(paste(cancer, "/countTable.RData", sep=""))
cds = newCountDataSet(round(countTable), tumor_info$type)
cds = estimateSizeFactors(cds)
cds = estimateDispersions(cds)

normalized_cds_all = counts( cds, normalized=TRUE )
tumor_info = subset(tumor_info, rownames(tumor_info) %in% colnames(normalized_cds_all))
#load(paste(cancer, "/cds.RData", sep=""))

load(paste("sigs", '/', cancer, '_sig_deseq_final.RData', sep='')) #tcga_deseq LIHC_sig_deseq_QC
res$GeneID = sapply((res$id), function(id){unlist(strsplit(id, "\\|"))[2]})
res$symbol = sapply((res$id), function(id){unlist(strsplit(id, "\\|"))[1]})
dz_signature = subset(res, !is.na(padj) & !is.na(id) & id !='?' & padj < 1E-3 & abs(log2FoldChange) > 2 & abs(log2FoldChange) != Inf )
dz_sig = subset(dz_signature, select=c("GeneID", "log2FoldChange"))

load("~/Documents/stanford/lincs/data/lincs_signatures_cmpd_landmark.RData")
lincs_sig_info = read.csv("~/Documents/stanford/lincs/data/lincs_sig_info.csv")

if (cancer == "LIHC"){
  cell_lines = c("HEPG2", "HUH7")
}else if (cancer == "BRCA"){
  cell_lines = c("MCF7", "BT20", "SKBR3")
}else if (cancer == "COAD"){
  cell_lines = c("HT29", "LOVO", "SNUC4", "SNUC5","SW948")
}
lincs_sig_info = subset(lincs_sig_info, cell_id %in% cell_lines)

#remove low expression hits
d = DGEList(counts=countTable)
d$samples = tumor_info
keep <- rowSums(cpm(d$counts)>1) >= as.numeric(sort(table(tumor_info$type))[1])
countTable = countTable[keep,]
gene_ids = sapply(rownames(normalized_cds_all), function(id){unlist(strsplit(id, "\\|"))[2]})

matched_patients = unique(tumor_info$bcr_patient_barcode[tumor_info$type == "tumor"])

#non tumor gene expression
patient_non_tumor = rownames(tumor_info[tumor_info$type == "non-tumor",])
non_tumor_sig = apply(normalized_cds_all[, patient_non_tumor], 1, mean)

patient_drug_cor = NULL
patient_sigs = data.frame(id = rownames(normalized_cds_all))
for (patient in matched_patients){
   patient_tumor = rownames(tumor_info[tumor_info$bcr_patient_barcode == patient & tumor_info$type == "tumor",])
   tumor_sig = normalized_cds_all[, patient_tumor[1]]
   
   #compare to non tumors
   patient_sig = log(tumor_sig+1) - log(non_tumor_sig + 1)
   patient_sigs = cbind(patient_sigs, patient_sig)
   
   #one patient may have multiple drugs
   drugs = c("niclosamide")
   
   drug_instances = subset(lincs_sig_info, tolower(lincs_sig_info$pert_iname) %in% tolower(as.character(drugs)) &
                             pert_type == "trt_cp" & is_gold ==1, select=c("id"))
   
   if (nrow(drug_instances)>0){
    drug_sigs = data.frame(lincs_signatures[,colnames(lincs_signatures) %in% c(as.character(drug_instances$id))])
    
    patient_drug = merge(data.frame(gene_id = gene_ids, patient_sig), data.frame(gene_id = rownames(drug_sigs), drug_sigs), by="gene_id")
    patient_drug = subset(patient_drug, gene_id %in% dz_sig$GeneID)
    cors = cor(patient_drug[,2:ncol(patient_drug)], method="spearman")
    cor = median(cors[1, 2:ncol(cors)]) #take the median of all correlations
   }else{
     cor = NA
   }
   patient_drug_cor = c(patient_drug_cor, cor)
}

patient_drug= data.frame(patient_drug_cor, matched_patients)
patient_drug = merge(patient_drug, tumor_info, by.x="matched_patients", by.y="bcr_patient_barcode")
patient_drug = subset(patient_drug, !is.na(patient_drug_cor) )
patient_drug = patient_drug[order(patient_drug$patient_drug_cor), ]

save(paste(cancer,"/niclosamide_patient.RData")
#median of niclosamide: -0.34; sorafenib: -0.22