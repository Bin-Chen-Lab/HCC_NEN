#scripts used to correct signficance between tumors and cell lines. Random samples were taken from expO
# the cutoff may be slightly different due to the random sampling.

library(GEOquery)

is_outlier <- function(cancer, cell_line_one_tumor_anno){
  #compute correlation between tumors and dz related cell line as well as other cell lines
  ccle_tcga = read.csv("ccle_cellline_tcga_mapping_updated.csv", stringsAsFactors=F)
  ccle_tcga_target_celllines = ccle_tcga$CCLE.name[ccle_tcga$tcga.tumor == cancer]
  cor_same = cell_line_one_tumor_anno$cor[cell_line_one_tumor_anno$CCLE.name  %in% ccle_tcga_target_celllines]
  cor_others = cell_line_one_tumor_anno$cor[!(cell_line_one_tumor_anno$CCLE.name  %in% ccle_tcga_target_celllines)]
  p = wilcox.test(cor_same, cor_others, alternative = 'greater')  
  return(p$p.value)
}


raw_matrix <- read.csv("raw/expo/GSE2109_series_matrix_cleaned.txt", sep="\t")
raw_matrix <- raw_matrix[, c(1,sample(2:ncol(raw_matrix), 1000))] #random take 1000 tumors

gpl <- getGEO(filename=("raw/expo/GPL570.soft"))
geneMappings <- Table(gpl)[,c("ENTREZ_GENE_ID","ID")]
names(geneMappings) <- c("GeneID","ID")
gpl <- geneMappings[!duplicated(geneMappings$ID),]

merge_matrix = merge(gpl, raw_matrix, by.x="ID", by.y="ID_REF")
merge_matrix = merge_matrix[,-c(1)]
merge_matrix = aggregate(. ~ GeneID, merge_matrix, mean)


load("ccle_meta_updated.RData")
load("ccle_gene_updated.RData")
iqr_gene = apply(ccle_gene[,-c(1)], 1, IQR)
gene_ids = ccle_gene$GeneID
varying_genes = gene_ids[(order(iqr_gene, decreasing=T))][1:num_varying_genes]
ccle_gene = subset(ccle_gene, GeneID %in% varying_genes)

cell_line_tumor = merge( ccle_gene, merge_matrix, by.x="GeneID", by.y= "GeneID")
cell_line_tumor_cor = cor(cell_line_tumor[,-c(1)], method="spearman")

tumors = colnames(merge_matrix)[-1]
celllines = colnames(ccle_gene)[-1]

cancers = cancer
cutoffs = NULL
for (cancer in cancers){
  tumor_outliers = NULL
  for (tumor in tumors){
    cell_line_one_tumor_cor = cell_line_tumor_cor[tumor, celllines ] 
    cell_line_one_tumor_cor = data.frame(sample = names(cell_line_one_tumor_cor), cor=cell_line_one_tumor_cor)
    cell_line_one_tumor_anno = merge(ccle_meta, cell_line_one_tumor_cor, by.y="sample", by.x="CCLE.name")
    cell_line_one_tumor_anno = subset(cell_line_one_tumor_anno, select=c("Cell.line.primary.name", "CCLE.name", "cor", "Histology","Hist.Subtype1" ,"Site.Primary"))  
    
    tumor_outliers = c(tumor_outliers, is_outlier(cancer, cell_line_one_tumor_anno))
  }
  cutoffs = c(cutoffs, sort(tumor_outliers)[50]) #the p value in the 50th out of 1000 (5%)
}

#cutoff used
cutoff = cutoffs[1]
cutoff