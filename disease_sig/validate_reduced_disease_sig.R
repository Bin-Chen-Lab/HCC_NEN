#validate reduced disease signatures using TCGA data
#create heatmap

library(pheatmap)
library(gplots)

landmark <- read.csv("raw/lincs/lincs_landmark.csv")

cancer <- "LIHC"
dz_padj_cutoff <- 1E-3
dz_fc_cutoff <- 2

load(paste(cancer,"/cds.RData", sep=""))
load(paste(cancer,"/tumor_info.RData", sep=""))

load(paste( cancer, "/", cancer, '_sig_deseq_final.RData', sep=''))
dz_signature <- subset(res, !is.na(padj) & !is.na(id) & id !='?' & padj < dz_padj_cutoff & abs(log2FoldChange) > dz_fc_cutoff & abs(log2FoldChange) != Inf )

normalized_cds_all <- counts( cds, normalized=TRUE )
normalized_cds_all <- log(normalized_cds_all + 1)

#remove outlier samples
normalized_cds_all <- normalized_cds_all[, colnames(normalized_cds_all) %in% rownames(tumor_info[tumor_info$type == "non-tumor" | (tumor_info$type == "tumor" & tumor_info$quality != "outlier"),])]

#GeneID <- sapply(rownames(normalized_cds_all), function(id){
#  as.character(unlist(strsplit(id, "\\|"))[2])
#})

GeneID <- sapply(dz_signature$id, function(id){
  as.character(unlist(strsplit(id, "\\|"))[2])
})

GeneID <- GeneID[GeneID %in% landmark$gene_id]

normalized_cds_sig <- normalized_cds_all[GeneID,]

annotation <- subset(tumor_info, select= c("tcga_barcode", "type"))
rownames(annotation) <- annotation$tcga_barcode
annotation$type <- as.factor(annotation$type)
annotation <- subset(annotation, select=c("type"))


Var1        <- c("lightblue", "green")
names(Var1) <- c("tumor", "non-tumor")
anno_colors <- list(type = Var1)

my.cols <- greenred(100) #brewer.pal(9, "Blues")
pheatmap(scale(normalized_cds_sig), col = my.cols, annotation = annotation, annotation_colors = anno_colors,
         show_colnames=F, legend=T, show_rownames=F, filename=paste(cancer, "/dz_sig/reduced_dz_sig_validation.pdf", sep="")
)


