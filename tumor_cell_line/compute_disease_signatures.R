# TCGA data analysis

library("RColorBrewer")
library("gplots")
library(DESeq)

###########
##process RNASEQ from TCGA
#quality control
#cross validation; validate with independent sets
cutoff = .05
data_from_gdac = T

meta_file = list.files(paste("mrna/",cancer,"/METADATA/UNC__IlluminaHiSeq_RNASeqV2//", sep=''), 'sdrf', full.names=T)[1]
tcga_meta = read.csv(meta_file, sep="\t", stringsAsFactors=F)

tcga_meta = subset(tcga_meta, Protocol.REF.4 == "unc.edu:RSEM_genes:IlluminaHiSeq_RNASeqV2:2")

patient_id <- sapply(tcga_meta$Comment..TCGA.Barcode., function(barcode){
  tags = unlist(strsplit(as.character(barcode), "-"))
  paste(tags[1:3], collapse="-")
})

sample_id <- sapply(tcga_meta$Comment..TCGA.Barcode., function(barcode){
  tags = unlist(strsplit(as.character(barcode), "-"))
  paste(tags[1:4], collapse="-")
})

tumor_type <- sapply(tcga_meta$Comment..TCGA.Barcode., function(barcode){
  tags = unlist(strsplit(as.character(barcode), "-"))
  id = as.numeric(paste(unlist(strsplit(tags[4], ""))[1:2], collapse=""))
  if (id <= 9){
    "tumor"
  }else if (id > 9 & id <= 19){
    "non-tumor"
  }else{
    "control"
  }
})


if (data_from_gdac){
  raw.data.header = read.csv(paste('firehose/gdac.broadinstitute.org_', cancer, '.Merge_rnaseqv2__illuminahiseq_rnaseqv2__unc_edu__Level_3__RSEM_genes__data.Level_3.2014051800.0.0/', cancer, '.rnaseqv2__illuminahiseq_rnaseqv2__unc_edu__Level_3__RSEM_genes__data.data.txt', sep=""), sep='\t', nrow=1, check.names=F)
  
  raw.data = read.csv(paste('firehose/gdac.broadinstitute.org_', cancer, '.Merge_rnaseqv2__illuminahiseq_rnaseqv2__unc_edu__Level_3__RSEM_genes__data.Level_3.2014051800.0.0/', cancer, '.rnaseqv2__illuminahiseq_rnaseqv2__unc_edu__Level_3__RSEM_genes__data.data.txt', sep=""), sep='\t', skip=1, check.names=F)
  genes = as.character(raw.data[,1])
  GeneID = sapply(genes, function(id){
    as.character(unlist(strsplit(id, "\\|"))[2])
  })
  
  
  countTable = raw.data[,colnames(raw.data) == 'raw_count']
  raw.data.header = colnames(raw.data.header)[colnames(raw.data) == 'raw_count']
  colnames(countTable) = raw.data.header 
  rownames(countTable) = genes #GeneID
  invalid_files = which(!as.character(tcga_meta$Comment..TCGA.Barcode.) %in% colnames(countTable))
  if (length(invalid_files) > 0){
    countTable = countTable[, as.character(tcga_meta$Comment..TCGA.Barcode.)[-c(invalid_files)]] #reorder sample
    patient_id = patient_id[-c(invalid_files)]
    tumor_type = tumor_type[-c(invalid_files)]
    sample_id = sample_id[-c(invalid_files)]
    tcga_barcode = tcga_meta$Comment..TCGA.Barcode.[-c(invalid_files)]
  }else{
    countTable = countTable[, as.character(tcga_meta$Comment..TCGA.Barcode.)] #reorder sample    
  }
}else{
  invalid_files = NULL
  countTable = data.frame()
  for (i in 1:nrow(tcga_meta)){
    source_id = tcga_meta$Derived.Data.File[i]
    source_path = paste("mrna/", cancer,"/RNASeqV2/UNC__IlluminaHiSeq_RNASeqV2/Level_3/", source_id, sep="")
    if (!file.exists(source_path)){
      print(source_path)
      invalid_files = c(invalid_files, i)
      next
    }
    raw_data = read.csv(source_path, sep="\t")
      
    if (i ==1 ){
      countTable = data.frame(name = raw_data[,2])    
      names(countTable) = as.character(tcga_meta$Comment..TCGA.Barcode.[i])
      rownames(countTable) = raw_data[,1]
    }else{
      countTable[,as.character(tcga_meta$Comment..TCGA.Barcode.[i])] = raw_data[,2]
    }
  }
  
  if (length(invalid_files)>0){
    patient_id = patient_id[-c(invalid_files)]
    tumor_type = tumor_type[-c(invalid_files)]
    sample_id = sample_id[-c(invalid_files)]
  }
}




#stratify patients
clinical = read.csv(paste("clinical/", cancer, "/Clinical/Biotab/nationwidechildrens.org_clinical_patient_",cancer,".txt", sep=""), sep="\t")
patient_phenotype = merge(data.frame(bcr_patient_barcode=patient_id, tcga_barcode), clinical, by.x=1, by.y="bcr_patient_barcode", all.x=T, sort=F)
rownames(patient_phenotype) = patient_phenotype$tcga_barcode
patient_phenotype = patient_phenotype[tcga_barcode,]

biospecimen = read.csv(paste("clinical/", cancer, "/Clinical/Biotab/nationwidechildrens.org_biospecimen_sample_", cancer, ".txt", sep=""), sep="\t")
sample_phenotype = merge(sample_id, biospecimen, by.x=1, by.y="bcr_sample_barcode", all.x=T, sort=F)


###
tumor_cell_all = read.csv(paste(cancer, "/tumor_cell_line/", comparison_gene_set, "/tumor_cell_all", cancer, ".csv", sep=""), stringsAsFactors=F)
tumor_cell_all_p = aggregate(outlier ~ patient_id + sample_id, tumor_cell_all, min)
tumor_cell_all_p_adj = p.adjust(tumor_cell_all_p$outlier, "fdr")
patient_outliers =  unique(tumor_cell_all_p$patient_id[tumor_cell_all_p_adj > cutoff])
sample_outliers =  unique(tumor_cell_all_p$sample_id[tumor_cell_all_p_adj > cutoff])

patient_phenotype$quality <- sapply(patient_id, function(id){
  ifelse(id %in% patient_outliers, T, F)
})

sample_phenotype$quality <- sapply(sample_id, function(id){
  ifelse(id %in% sample_outliers, T, F)
})


#############
#remove the outliers

#####
#compute DE genes
#remove unqaulified samples

cds = newCountDataSet( round(countTable[, !(sample_id %in% sample_outliers)], 0), tumor_type[!(sample_id %in% sample_outliers)] ) #only integer is accepted
cds = estimateSizeFactors( cds )
sizeFactors( cds )
head( counts( cds, normalized=TRUE ) )

cds = estimateDispersions( cds) #, method ="pooled", sharingMode="gene-est-only")
save(cds, file=paste(cancer,"/cds.RData", sep=""))

str( fitInfo(cds) )
plotDispEsts( cds )
head( fData(cds) )
res = nbinomTest( cds,  "non-tumor", "tumor" )
head(res)

res_subset = subset(res, abs(log2FoldChange) != Inf)
plotMA(res_subset , col = ifelse(res_subset$padj < 1E-10 & abs(res_subset$log2FoldChange) > 2, "red", "gray"))

symbol= sapply(res$id, function(id){
  as.character(unlist(strsplit(id, "\\|"))[1])
})
GeneID = sapply(res$id, function(id){
  as.character(unlist(strsplit(id, "\\|"))[2])
})

res = cbind(res, symbol, GeneID)

save(res, file= paste(cancer, "/", cancer, "tcga_deseq.RData", sep=""))



