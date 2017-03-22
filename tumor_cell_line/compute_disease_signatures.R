# TCGA data analysis

library("RColorBrewer")
library("gplots")
library(DESeq)

###########
load(paste(cancer,"/countTable.RData", sep=""))

#remove genes with weak signals
load(paste(cancer,"/tumor_info.RData", sep=""))
countTable = round(countTable)
#remove outliers
tumor_info = tumor_info[tumor_info$quality != "outlier",]
keep <- rowSums(countTable>1) >= as.numeric(sort(table(tumor_info$type))[1])

save(countTable, file=paste(cancer,"/countTable.RData", sep=""))

#####
#compute DE genes
#remove inconsistent samples
cds = newCountDataSet( round(countTable[keep, !(sample_id %in% sample_outliers)], 0),  tumor_type[!(sample_id %in% sample_outliers)] ) #only integer is accepted
cds = estimateSizeFactors( cds )
sizeFactors( cds )
head( counts( cds, normalized=TRUE ) )

#visualize pca
cdsB = estimateDispersions(cds, method = "blind")
vsd = varianceStabilizingTransformation(cdsB)
p = plotPCA(vsd, "condition")

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

save(res, file= paste(cancer, "/", cancer, "tcga_deseq_final.RData", sep=""))



