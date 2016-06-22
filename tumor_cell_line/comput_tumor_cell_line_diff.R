#compute differentially expressed genes between tumors and cell lines

library(DESeq)
library(pheatmap)
library("RColorBrewer")
layout(matrix(1))

load(paste(cancer,"/cds.RData", sep=""))
load(paste(cancer,"/tumor_info.RData", sep=""))
load("ccle_gene_updated.RData")
load("ccle_meta_updated.RData")

#raw count table
normalized_cds_all = counts( cds, normalized=TRUE )
tumor_samples = rownames(tumor_info[tumor_info$type == "tumor",])
non_tumor_samples = rownames(tumor_info[tumor_info$type == "non-tumor",])
tumor_normalized_cds = normalized_cds_all[, colnames(normalized_cds_all) %in% tumor_samples]

gene_ids = sapply(rownames(tumor_normalized_cds), function(gene){
  as.numeric(unlist(strsplit(gene, "\\|"))[2])
})
symbols = sapply(rownames(tumor_normalized_cds), function(gene){
  as.character(unlist(strsplit(gene, "\\|"))[1])
})

#ranked expression
tumor_normalized_cds = tumor_normalized_cds[gene_ids %in% ccle_gene$GeneID,]
for(i in 1:ncol(tumor_normalized_cds)){
  tumor_normalized_cds[,i] = rank(tumor_normalized_cds[,i])
}

#find hcc related cell lines
tumor_cell_mapping = read.csv("ccle_cellline_tcga_mapping_updated.csv")
tumor_cell_mapping = subset(tumor_cell_mapping, select=c("CCLE.name", "tcga.tumor", "tumor_type_name", "Site.Primary"))
names(tumor_cell_mapping) = c("CCLE.name", "tumor_type", "cell_line_type_name", "primary_site")
tumor_cell_mapping_subset = subset(tumor_cell_mapping, tumor_type == cancer)
ccle_meta_hcc = subset(ccle_meta, ccle_meta$CCLE.name %in% tumor_cell_mapping_subset$CCLE.name)

ccle_gene_subset = ccle_gene[ccle_gene$GeneID %in% gene_ids, colnames(ccle_gene) %in% c("GeneID", as.character(ccle_meta_hcc$CCLE.name))]

#rank expression for cell lines
for(i in 2:ncol(ccle_gene_subset)){
  ccle_gene_subset[,i] = rank(ccle_gene_subset[,i])
}

p_value_t_test = NULL
p_value_wilcox = NULL
fc = NULL
gene_ids = NULL
symbols = NULL
for (enzyme in rownames(tumor_normalized_cds)){
  gene_id = as.numeric(unlist(strsplit(enzyme, "\\|"))[2])
  gene_ids = c(gene_ids, gene_id)
  symbols = c(symbols, unlist(strsplit(enzyme, "\\|"))[1])
  tumor = tumor_normalized_cds[enzyme,]
  cell = ccle_gene_subset[ccle_gene_subset$GeneID == gene_id, -c(1)]
  if (nrow(cell)>0){
    p_value_t_test = c(p_value_t_test, t.test(tumor, as.numeric(cell))$p.value)
    p_value_wilcox = c(p_value_wilcox, wilcox.test(tumor, as.numeric(cell))$p.value)
    fc = c(fc, mean(tumor)/mean(as.numeric(cell)))
  }else{
    p_value_t_test = c(p_value_t_test, 1)
    p_value_wilcox = c(p_value_wilcox, 1)
    fc = c(fc, 0)
  }
}

enzyme_p = data.frame(rownames(tumor_normalized_cds), p_value_t_test,p_value_wilcox, fc , symbols, gene_ids)
enzyme_p$q = p.adjust(enzyme_p$p_value_t_test, "fdr")

enzyme_p = enzyme_p[order(enzyme_p$fc),]
write.csv(enzyme_p, paste(cancer, "/tumor_cell_line/genes_tumor_cell_line.csv", sep=""))

DE_genes = enzyme_p[(enzyme_p$fc>4 | enzyme_p$fc<0.25) & enzyme_p$q < 0.01,  ]

write.csv(DE_genes, paste(cancer, "/tumor_cell_line/DE_genes_tumor_cell_line.csv", sep=""))



#visualize differentially expressed genes in tumors or cell lines
gene_ids = ccle_gene$GeneID
ccle_gene_scaled = as.data.frame((ccle_gene[,-c(1)]))
ccle_gene_scaled$GeneID = gene_ids

tumor_normalized_cds_log = log(tumor_normalized_cds + 1)

my.col = colorRampPalette(brewer.pal(8,"YlOrRd"))(100)

#visualize all differentially expressed genes in TCGA
gene_ids = sapply(rownames(tumor_normalized_cds_log), function(gene){
  as.numeric(unlist(strsplit(gene, "\\|"))[2])
})


DE_genes = enzyme_p$gene_ids[(enzyme_p$fc>4 | enzyme_p$fc<0.25) & enzyme_p$q < 0.01 ]
tumor_normalized_cds_log_subset = tumor_normalized_cds_log[gene_ids %in% DE_genes, ]
symbols = sapply(rownames(tumor_normalized_cds_log_subset), function(gene){
  as.character(unlist(strsplit(gene, "\\|"))[1])
})
rownames(tumor_normalized_cds_log_subset) = symbols
pheatmap((tumor_normalized_cds_log_subset), color= my.col, cellheight=8, cellwidth=10, cex=0.8,file=paste(cancer, "/tumor_cell_line/DE_tumor_expression.pdf", sep=""))

#visualize all differentially expressed genes
DE_genes = enzyme_p$gene_ids[(enzyme_p$fc>4 | enzyme_p$fc<0.25) & enzyme_p$q < 0.01 ]

ccle_gene_subset = ccle_gene_scaled[ccle_gene_scaled$GeneID %in% DE_genes, colnames(ccle_gene_scaled) %in% c("GeneID", as.character(ccle_meta_hcc$CCLE.name))]
gene_ids = ccle_gene_subset$GeneID
ccle_gene_subset = ccle_gene_subset[, !(colnames(ccle_gene_subset) %in% c("GeneID"))]
rownames(ccle_gene_subset) = merge(gene_ids, enzyme_p, by.x=1, by.y="gene_ids", sort=F)$symbols
pheatmap((ccle_gene_subset), color= my.col, cellheight=8, cellwidth=8, cex=0.8,file=paste(cancer, "/tumor_cell_line/DE_cellline_expression.pdf", sep=""))



#different expressed enzymes
load(paste(cancer, "/", cancer, "tcga_deseq.RData", sep=""))

#overlap DE genes between tumor & non tumor and tumor & cell lines
phyper.test <- function(list1, list2, background){
  #gene list1: balls drawn, gene_list2: white balls in background
  
  #first step should remove the genes not shown in background, the step is moved to the parent step to increase efficiency
  list1 <- unique(list1[list1 %in% background])
  list2 <- unique(list2[list2 %in% background])
  
  q = sum(list1 %in% list2)
  m = length(list2)
  n = length(background) - length(list2)
  k = length(list1)
  # p = 1 - phyper(q,m,n,k)
  p = phyper(q,m,n,k,lower.tail=F,log.p=T)
  
  return(exp(p))
}

tumor_cell_line = read.csv(paste(cancer, "/tumor_cell_line/genes_tumor_cell_line.csv", sep=""))
DE_tumor_cell_line = subset(tumor_cell_line, (fc < 0.25 ) & q < 0.01)
dim(DE_tumor_cell_line)
DE_tumor_non_tumor = subset(res, !is.na(padj) & !is.na(id) & id !='?' & padj < 0.01 & abs(log2FoldChange) > 2 & abs(log2FoldChange) != Inf )
dim(DE_tumor_non_tumor)

common_genes = merge(DE_tumor_non_tumor, DE_tumor_cell_line, by.x="GeneID", by.y="gene_ids")
dim(common_genes)
background_genes = intersect(as.numeric(tumor_cell_line$gene_ids), as.numeric(res$GeneID))                 
phyper.test(as.numeric(DE_tumor_non_tumor$GeneID), as.numeric(DE_tumor_cell_line$gene_ids), background_genes )

write.csv(DE_tumor_non_tumor, paste(cancer, "/tumor_cell_line/genes_tumor_non_tumor.csv", sep=""))
write.csv(subset(common_genes, select=c("GeneID", "symbol")), paste(cancer, "/tumor_cell_line/common_down_tumors_DE.csv", sep=""))



