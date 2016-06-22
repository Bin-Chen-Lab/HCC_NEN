#compare expressions between TCGA tumors  and cell lines; detect the suspious samples; visualize the samples

library("RColorBrewer")
library("gplots")
library("ggplot2")
library(DESeq)
library("RMySQL")
library("ROCR")
library("beanplot")
mysql_drvr <-dbDriver("MySQL")
#con <- dbConnect(mysql_drvr,group="client",host="buttelab-db1.stanford.edu",dbname="proj_lincs")


########
#functions

is_outlier <- function(cancer, cell_line_one_tumor_anno){
  #compute correlation between tumors and dz related cell line as well as other cell lines
  ccle_tcga = read.csv("ccle_cellline_tcga_mapping_updated.csv", stringsAsFactors=F)
  ccle_tcga_target_celllines = ccle_tcga$CCLE.name [ccle_tcga$tcga.tumor == cancer]
  cor_same = cell_line_one_tumor_anno$cor[cell_line_one_tumor_anno$CCLE.name  %in% ccle_tcga_target_celllines]
  cor_others = cell_line_one_tumor_anno$cor[!(cell_line_one_tumor_anno$CCLE.name  %in% ccle_tcga_target_celllines)]
  p = wilcox.test(cor_same, cor_others, alternative = 'greater')  
  return(p$p.value)
}

is_outlier3 <- function(cancer, cell_line_one_tumor_anno){
  ccle_tcga = read.csv("ccle_cellline_tcga_mapping_updated.csv", stringsAsFactors=F)
  ccle_tcga_target_celllines = ccle_tcga$CCLE.name [ccle_tcga$tcga.tumor == cancer]
  cor_same = cell_line_one_tumor_anno$cor[cell_line_one_tumor_anno$CCLE.name  %in% ccle_tcga_target_celllines]
  cor_others = cell_line_one_tumor_anno$cor[!(cell_line_one_tumor_anno$CCLE.name  %in% ccle_tcga_target_celllines)]
  pred <- prediction(c(cor_same, cor_others), c(rep(T, length(cor_same)), rep(F, length(cor_others))))
  perf <- performance(pred, measure = "auc") 
  
  return(perf@y.values[[1]])
}

bean.plot <- function(tumor_cell_cor, cancer, type){
  group = sapply(1:nrow(tumor_cell_cor), function(id){
    if (tumor_cell_cor[id, "histology_class"] == "others"){
      paste(tumor_cell_cor[id, "sample_id"], 1)
    }else{
      paste(tumor_cell_cor[id, "sample_id"], 2)
    }
  })
  
  tumor_cell_cor$group = group
  
  pdf(paste(cancer, "/tumor_cell_line/", comparison_gene_set, "/", cancer, "_", type, "_beanplot.pdf",sep=""))
  par(mar = c(10, 4, 4, 2) + 0.1, las=2)
  beanplot(cor ~ group, data = tumor_cell_cor, ll = 0.15,
           main = "", ylab = "correlation", side = "both",
           border = NA, col = list(c("grey","green","green"), c("grey", "red","red")), 
           beanlines = "median", overalllin="median", method="jitter")
  legend("bottomleft", fill = c("green", "red"),
         legend = c("others", paste(cancer_name,"related")))
  dev.off()
}

box.plot <- function(tumor_cell_cor, cancer, type){
  pdf(paste(cancer, "/tumor_cell_line/", comparison_gene_set, "/", cancer, "_", type, "boxplot.pdf",sep=""))
  #random select 5 samples
  p <- ggplot(tumor_cell_cor, aes(factor(sample_id), cor, color= histology_class))
  print(p +   geom_boxplot(outlier.colour = "grey", notch=F, outlier.shape = NA) + geom_jitter() +   
    ylab("correlation") +
    xlab("sample") +  theme(axis.text.x = element_text(angle = 45, hjust = 1)))
  dev.off()
}

plot.sample <- function(tumor_cell_all_subset, cancer, sample1, type){
  dir.create(paste(cancer, "/tumor_cell_line/", comparison_gene_set, "/", "samples_", cancer,  "_", type, sep=""), showWarnings = FALSE)
  tumor_cell_all_outliers_sample = subset(tumor_cell_all_subset, as.character(sample_id) == sample1)
  top_tumor_types = as.character(get.top.cor.tumor.type(tumor_cell_all_outliers_sample, sample1))
  
  tumor_cell_cor= subset(tumor_cell_all_outliers, tumor_type %in% unique(c(cancer, top_tumor_types)))
  
  #order based on median cor
  tumor_cell_cor_merged = aggregate(cor ~ tumor_type_name, tumor_cell_cor, median)
  tumor_cell_cor_merged = tumor_cell_cor_merged[order(tumor_cell_cor_merged$cor, decreasing=T), ]
  tumor_cell_cor$tumor_type_name = factor(tumor_cell_cor$tumor_type_name, levels = tumor_cell_cor_merged$tumor_type_name)
  
  pdf(paste(cancer, "/tumor_cell_line/", comparison_gene_set, "/", "samples_", cancer, "_", type, "/", sample1, ".pdf",sep=""))
  
  p <- ggplot(tumor_cell_cor, aes(tumor_type_name, cor))
  print(p +   geom_boxplot(outlier.colour = "grey", notch=F, outlier.shape = NA) + geom_jitter() +   
          ylab("correlation") +
          xlab("") +  theme(axis.text.x = element_text(angle = 45, hjust = 1)))
  dev.off()
}

get.top.cor.tumor.type <- function(tumor_cell_all_outliers_sample, sample1){
  tumor_types = unique(tumor_cell_all_outliers_sample$tumor_type)
  
  p_values = NULL
  for (tumor_type in tumor_types){
    
    ccle_tcga_target_celllines = tumor_cell_all_outliers_sample$CCLE.name[tumor_cell_all_outliers_sample$tumor_type == tumor_type]
    cor_same = tumor_cell_all_outliers_sample$cor[tumor_cell_all_outliers_sample$CCLE.name %in% ccle_tcga_target_celllines]
    cor_others = tumor_cell_all_outliers_sample$cor[!(tumor_cell_all_outliers_sample$CCLE.name %in% ccle_tcga_target_celllines)]
    p = wilcox.test(cor_same, cor_others, alternative = 'greater')  
    p_values = c(p_values, p$p.value)
  }
  
  tumor_type_p = data.frame(tumor_type = tumor_types, p = p_values)
  tumor_type_p = tumor_type_p[order(tumor_type_p$p),]
  return (tumor_type_p$tumor_type[1:5])
}


###########
##process RNASEQ from TCGA
#quality control
#cross validation; validate with independent sets


if (data_from_gdac){
  raw.data.header = read.csv(paste('firehose/gdac.broadinstitute.org_', cancer, '.Merge_rnaseqv2__illuminahiseq_rnaseqv2__unc_edu__Level_3__RSEM_genes__data.Level_3.2014051800.0.0/', cancer, '.rnaseqv2__illuminahiseq_rnaseqv2__unc_edu__Level_3__RSEM_genes__data.data.txt', sep=""), sep='\t', nrow=1, check.names=F)
  
  raw.data = read.csv(paste('firehose/gdac.broadinstitute.org_', cancer, '.Merge_rnaseqv2__illuminahiseq_rnaseqv2__unc_edu__Level_3__RSEM_genes__data.Level_3.2014051800.0.0/', cancer, '.rnaseqv2__illuminahiseq_rnaseqv2__unc_edu__Level_3__RSEM_genes__data.data.txt', sep=""), sep='\t', skip=1, check.names=F)
  genes = as.character(raw.data[,1])
  GeneID = sapply(genes, function(id){
    as.character(unlist(strsplit(id, "\\|"))[2])
  })
  
  
  countTable = raw.data[,colnames(raw.data) == 'scaled_estimate']
  raw.data.header = colnames(raw.data.header)[colnames(raw.data) == 'scaled_estimate']
  colnames(countTable) = raw.data.header 
  rownames(countTable) = GeneID
}else{
  meta_file = paste("mrna/",cancer,"/METADATA/UNC__IlluminaHiSeq_RNASeqV2//unc.edu_", cancer, ".IlluminaHiSeq_RNASeqV2.1.10.0.sdrf.txt", sep='')
  if (!file.exists(meta_file)){
    meta_file = paste("mrna/",cancer,"/METADATA/UNC__IlluminaHiSeq_RNASeqV2//unc.edu_", cancer, ".IlluminaHiSeq_RNASeqV2.1.11.0.sdrf.txt", sep='')
    if (!file.exists(meta_file)){
      print(paste(meta_file, "does not exists"))
    }
  }
  
  tcga_meta = read.csv(meta_file, sep="\t")
  tcga_meta = subset(tcga_meta, Protocol.REF.4 == "unc.edu:RSEM_genes:IlluminaHiSeq_RNASeqV2:2")
  
  #get expression
  countTable = data.frame()
  invalid_files = NULL
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
      countTable = data.frame(name = raw_data[,3])    #use scale estimate
      names(countTable) = as.character(tcga_meta$Comment..TCGA.Barcode.[i])
      rownames(countTable) = raw_data[,1] #gene id
    }else{
      countTable[,as.character(tcga_meta$Comment..TCGA.Barcode.[i])] = raw_data[,3]
    }
  }
  
  
  GeneID = sapply(rownames(countTable), function(id){
    as.character(unlist(strsplit(id, "\\|"))[2])
  })
}

if (!file.exists(cancer)){
  dir.create(cancer)
}
if (!file.exists(paste(cancer, "/tumor_cell_line", sep=""))){
  dir.create(paste(cancer, "/tumor_cell_line", sep=""))
}

##compare with expression in ccle
load("ccle_meta_updated.RData")
load("ccle_gene_updated.RData")

if (comparison_gene_set == "sigs"){
  load(paste('sigs/', cancer, '_sig_deseq_final.RData', sep=''))
  res$GeneID = sapply(res$id, function(id){
    as.character(unlist(strsplit(id, "\\|"))[2])
  })
  dz_signature = subset(res, !is.na(padj) & !is.na(id) & id !='?' & padj < 1E-12 & abs(log2FoldChange) > 2 & abs(log2FoldChange) != Inf )
  selected_genes = dz_signature$GeneID
}else if (comparison_gene_set == "varying5k"){
  iqr_gene = apply(ccle_gene[,-c(1)], 1, IQR)
  gene_ids = ccle_gene$GeneID
  varying_genes = gene_ids[(order(iqr_gene, decreasing=T))][1:num_varying_genes] 
  selected_genes = varying_genes
}else if (comparison_gene_set == "metabolism"){
  enzymes = read.csv(paste(cancer, "/tumor_cell_line/enyzmes.csv", sep=""))
  selected_genes = enzymes$GeneID
}else if (comparison_gene_set == "meta_sigs"){
  dz_signature = read.csv("../code/meta_analysis2/dz_sig_meta_analysis.txt", sep="\t")
  selected_genes = dz_signature$GeneID
}else if (comparison_gene_set == "varying5k_tumor"){
  iqr_gene = apply(countTable, 1, IQR)
  gene_ids = as.numeric(rownames(countTable))
  varying_genes = gene_ids[(order(iqr_gene, decreasing=T))][1:num_varying_genes] 
  selected_genes = varying_genes  
}else{
  stop("error: need comparison gene set")
}

ccle_gene = subset(ccle_gene, GeneID %in% selected_genes)    

tumor_sig = data.frame(GeneID, countTable)
cell_line_tumor = merge( ccle_gene, tumor_sig, by="GeneID")
cell_line_tumor_cor = cor(cell_line_tumor[,-c(1)], method="spearman")

tumors = colnames(tumor_sig)[-1]
celllines = colnames(ccle_gene)[-1]

tumor_outliers = NULL
tumor_cell_all = data.frame()
for (tumor in tumors){
  cell_line_one_tumor_cor = cell_line_tumor_cor[tumor, celllines ] 
  cell_line_one_tumor_cor = data.frame(sample = names(cell_line_one_tumor_cor), cor=cell_line_one_tumor_cor)
  cell_line_one_tumor_anno = merge(ccle_meta, cell_line_one_tumor_cor, by.y="sample", by.x="CCLE.name")
  cell_line_one_tumor_anno = subset(cell_line_one_tumor_anno, select=c("Cell.line.primary.name", "CCLE.name", "cor", "Histology","Hist.Subtype1" ,"Site.Primary"))  
  cell_line_one_tumor_anno$barcode = paste(unlist(strsplit(tumor, '\\.')), collapse="-")
  tags = unlist(strsplit(as.character(tumor), "\\."))
  cell_line_one_tumor_anno$patient_id =   paste(tags[1:3], collapse="-")
  cell_line_one_tumor_anno$sample_id =   paste(tags[1:4], collapse="-")
  #tumor_outliers = c(tumor_outliers, is_outlier(cancer, cell_line_one_tumor_anno))
  cell_line_one_tumor_anno$outlier = is_outlier(cancer, cell_line_one_tumor_anno)
  tumor_cell_all = rbind(tumor_cell_all, cell_line_one_tumor_anno)
}

tumor_cell_all = tumor_cell_all[order(tumor_cell_all$outlier,tumor_cell_all$barcode, tumor_cell_all$cor, decreasing=T),]

if (!file.exists(paste(cancer, "/tumor_cell_line/", comparison_gene_set, "/", sep=""))){
  dir.create(paste(cancer, "/tumor_cell_line/", comparison_gene_set, "/", sep=""))
}
write.csv(tumor_cell_all, paste(cancer, "/tumor_cell_line/", comparison_gene_set, "/", "tumor_cell_all", cancer, ".csv", sep=""))

###
#detect outliers
tumor_cell_all_p = aggregate(outlier ~ patient_id + sample_id, tumor_cell_all, min)
tumor_cell_all_p_adj = p.adjust(tumor_cell_all_p$outlier, "fdr")
patient_outliers =  unique(tumor_cell_all_p$patient_id[tumor_cell_all_p_adj > cutoff])
sample_outliers =  unique(tumor_cell_all_p$sample_id[tumor_cell_all_p_adj > cutoff])

write(patient_outliers, paste(cancer,"/tumor_cell_line/", comparison_gene_set, "/", "outlier_", cancer, "_outlier.txt", sep=""), ncolumn=1)

#visualize outliers
ccle_tcga = read.csv("ccle_cellline_tcga_mapping_updated.csv")
ccle_tcga_target_celllines = ccle_tcga$CCLE.name[ccle_tcga$tcga.tumor == cancer]

histology_class <- sapply(tumor_cell_all$CCLE.name, function(cell_line){
  if (cell_line %in% ccle_tcga_target_celllines){
    paste("HCC", "related")
  }else{
    "others"
  }
})

tumor_cell_all$histology_class <-  histology_class

tumor_cell_all_outliers = tumor_cell_all[ as.character(tumor_cell_all$sample_id) %in% sample(as.character(sample_outliers), min(10, length(sample_outliers))),]
bean.plot(tumor_cell_all_outliers, cancer, "bad")
box.plot(tumor_cell_all_outliers, cancer, "bad")

good_samples = sample(unique(tumor_cell_all$sample_id[tumor_cell_all$outlier < cutoff]), min(10, length(unique(tumor_cell_all$sample_id[tumor_cell_all$outlier < cutoff]))))
tumor_cell_all_good = tumor_cell_all[ as.character(tumor_cell_all$sample_id) %in% as.character(good_samples),]
bean.plot(tumor_cell_all_good, cancer, "good")
box.plot(tumor_cell_all_good, cancer, "good")

#find top 5 enriched cell type
tumor_cell_mapping = read.csv("ccle_cellline_tcga_mapping_updated.csv")
tumor_cell_mapping = subset(tumor_cell_mapping, select=c("CCLE.name", "tcga.tumor", "tumor_type_name"))
names(tumor_cell_mapping) = c("CCLE.name", "tumor_type", "tumor_type_name")

tumor_cell_all_outliers = merge(tumor_cell_all_outliers, tumor_cell_mapping,  by.x="CCLE.name")
samples = unique(tumor_cell_all_outliers$sample_id)
for(sample1 in samples){
  plot.sample(tumor_cell_all_outliers, cancer, sample1, "bad")
}

tumor_cell_all_good = merge(tumor_cell_all_good, tumor_cell_mapping,  by.x="CCLE.name")
samples = unique(tumor_cell_all_good$sample_id)
for(sample1 in samples){
  plot.sample(tumor_cell_all_good, cancer, sample1, "good")
}

############
#clinical features
#use raw count instead of scaled estimate, 
raw.data.header = read.csv(paste('firehose/gdac.broadinstitute.org_', cancer, '.Merge_rnaseqv2__illuminahiseq_rnaseqv2__unc_edu__Level_3__RSEM_genes__data.Level_3.2014051800.0.0/', cancer, '.rnaseqv2__illuminahiseq_rnaseqv2__unc_edu__Level_3__RSEM_genes__data.data.txt', sep=""), sep='\t', nrow=1, check.names=F)

countTable = raw.data[,colnames(raw.data) == 'raw_count']
raw.data.header = colnames(raw.data.header)[colnames(raw.data) == 'raw_count']
colnames(countTable) = raw.data.header 
rownames(countTable) = GeneID

meta_file = list.files(paste("mrna/",cancer,"/METADATA/UNC__IlluminaHiSeq_RNASeqV2//", sep=''), 'sdrf', full.names=T)[1]
tcga_meta = read.csv(meta_file, sep="\t", stringsAsFactors=F)

tcga_meta = subset(tcga_meta, Protocol.REF.4 == "unc.edu:RSEM_genes:IlluminaHiSeq_RNASeqV2:2" | Protocol.REF.4 == "unc.edu:RSEM_genes:IlluminaHiSeq_RNASeqV2:3")

tcga_barcode = tcga_meta$Comment..TCGA.Barcode.

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

invalid_files = which(!as.character(tcga_meta$Comment..TCGA.Barcode.) %in% colnames(countTable))
if (length(invalid_files) > 0){
  countTable = countTable[, as.character(tcga_meta$Comment..TCGA.Barcode.)[-c(invalid_files)]] #reorder sample
  patient_id = patient_id[-c(invalid_files)]
  tumor_type = tumor_type[-c(invalid_files)]
  sample_id = sample_id[-c(invalid_files)]
  tcga_barcode = tcga_barcode[-c(invalid_files)]
}else{
  countTable = countTable[, as.character(tcga_meta$Comment..TCGA.Barcode.)] #reorder sample    
}

#tricky here... when merge, if the merged column in x is duplicated, the order in the new merged data frame will be changed
clinical = read.csv(paste("clinical/", cancer, "/Clinical/Biotab/nationwidechildrens.org_clinical_patient_",cancer,".txt", sep=""), sep="\t", stringsAsFactors=F)
patient_phenotype = merge(data.frame(bcr_patient_barcode=patient_id, tcga_barcode), clinical, by.x=1, by.y="bcr_patient_barcode", all.x=T, sort=F)
rownames(patient_phenotype) = patient_phenotype$tcga_barcode
patient_phenotype = patient_phenotype[tcga_barcode,]

biospecimen = read.csv(paste("clinical/", cancer, "/Clinical/Biotab/nationwidechildrens.org_biospecimen_sample_", cancer, ".txt", sep=""), sep="\t")
sample_phenotype = merge(sample_id, biospecimen, by.x=1, by.y="bcr_sample_barcode", all.x=T, sort=F)

patient_phenotype$quality <- sapply(patient_id, function(id){
  ifelse(id %in% patient_outliers, "outlier", "")
})

sample_phenotype$quality <- sapply(sample_id, function(id){
  ifelse(id %in% sample_outliers, "outlier", "")
})

tumor_info = data.frame(type=tumor_type, sample_phenotype, patient_phenotype)
save(tumor_info, file=paste(cancer,"/tumor_info.RData", sep=""))

cdsFull = newCountDataSet( round(countTable, 0), data.frame(type=tumor_type, sample_phenotype) ) #only integer is accepted
cdsFull = estimateSizeFactors( cdsFull )
cdsFullBlind = estimateDispersions( cdsFull, method = "blind" )
vsdFull = varianceStabilizingTransformation( cdsFullBlind )

pdf(paste(cancer, "/tumor_cell_line/", comparison_gene_set, "/", cancer, "pca_outlier.pdf", sep=""))
  print(plotPCA(vsdFull, intgroup=c( "type", "quality"), ntop = num_varying_genes))
dev.off()

pdf(paste(cancer, "/tumor_cell_line/", comparison_gene_set, "/",cancer, "pca.pdf", sep=""))
  print(plotPCA(vsdFull, intgroup=c( "type"), ntop = num_varying_genes))
dev.off()

cdsFullcounts = counts(cdsFullBlind, normalized = T)
cdsFullcounts = log(cdsFullcounts + 1)

#cluster visualization
my.cols = sapply(colnames(cdsFullcounts), function(sample_name){
  
  if (tumor_info[sample_name, "type"] == "tumor" &  tumor_info[sample_name, "quality"] == "outlier" ){
    1
  }else if (tumor_info[sample_name, "type"] == "tumor" &  tumor_info[sample_name, "quality"] != "outlier" ){
    2
  }else if (tumor_info[sample_name, "type"] == "non-tumor" &  tumor_info[sample_name, "quality"] == "outlier" ){
    3
  }else if (tumor_info[sample_name, "type"] == "non-tumor" &  tumor_info[sample_name, "quality"] != "outlier" ){
    4
  }else{
    5
  }
})

library(sparcl)
  hc <- hclust(dist(t(cdsFullcounts)),method="complete")
  pdf(paste(cancer, "/tumor_cell_line/",comparison_gene_set, "/",cancer, "_TCGA_cluster.pdf", sep=""))
    ColorDendrogram(hc, y= my.cols,main="",branchlength=25)
dev.off()


