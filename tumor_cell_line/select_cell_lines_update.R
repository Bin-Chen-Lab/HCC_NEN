#analyze tumor vs cell line correlations, pick good cell lines

library(pheatmap)
library(sparcl)
library(gplots)
library(ggplot2)
require(Heatplus)
library(RColorBrewer)
draw_colnames_45 <- function (coln, ...) {
  library(grid)
  m = length(coln)
  x = (1:m)/m - 1/2/m
  grid.text(coln, x = x, y = unit(0.96, "npc"), vjust = .5, 
            hjust = 1, rot = 45, gp = gpar(...)) ## Was 'hjust=0' and 'rot=270'
}


load(paste(cancer,"/tumor_info.RData", sep=""))
tumor_cell_mapping = read.csv("ccle_cellline_tcga_mapping_updated.csv")
tumor_cell_mapping = subset(tumor_cell_mapping, select=c("CCLE.name", "tcga.tumor", "tumor_type_name", "Site.Primary"))
names(tumor_cell_mapping) = c("CCLE.name", "tumor_type", "cell_line_type_name", "primary_site")
tumor_samples = tumor_info$tcga_barcode[tumor_info$type == "tumor" & tumor_info$quality != "outlier"]
  
tumor_cell_all = read.csv(paste(cancer, "/tumor_cell_line/", comparison_gene_set, "/tumor_cell_all", cancer, ".csv", sep=""), stringsAsFactors=F)
tumor_cell_all = merge(tumor_cell_all, tumor_cell_mapping,  by="CCLE.name")
tumor_cell_all = subset(tumor_cell_all, barcode %in% as.character(tumor_samples))

#cell line order
cor_cutoff = qnorm(0.05, mean(tumor_cell_all$cor), sd(tumor_cell_all$cor), lower.tail = F)

#tumor cell line network
tumor_cell_network = subset(tumor_cell_all, tumor_type == cancer & cor > cor_cutoff)
write.table(tumor_cell_network, paste(cancer, "/tumor_cell_line/", comparison_gene_set, "/tumor_cell_line_network.txt", sep=""), row.names=F, col.names=T, quote=F, sep="\t")

cell_line_order = aggregate(cor ~ CCLE.name, tumor_cell_all, median)
cell_line_order = cell_line_order[order(cell_line_order$cor, decreasing=T),]
head(cell_line_order, 10)

qnorm(0.01, mean(cell_line_order$cor), sd(cell_line_order$cor), lower.tail = F)

#sample order
sample_order = aggregate(cor ~ sample_id, tumor_cell_all, median)
sample_order = sample_order[order(sample_order$cor, decreasing=T),]
head(sample_order, 10)
tumor_info_cor = merge(tumor_info, sample_order, by.x= "x", by.y="sample_id")
tumor_info_cor = tumor_info_cor[order(tumor_info_cor$cor, decreasing=T),]
p <- ggplot(tumor_info_cor, aes(race, cor)) 
p + geom_boxplot()
#clinical_stage tumor_grade

#overall relations between tumors and cell lines
cell_lines = unique(tumor_cell_all$CCLE.name)
samples = unique(tumor_cell_all$sample_id)
cell_tumor_matrix = matrix(NA, nrow=length(cell_lines), ncol=length(samples), dimnames=list(cell_lines, samples))
for(i in 1:nrow(tumor_cell_all)){
  cell_tumor_matrix[tumor_cell_all$CCLE.name[i],tumor_cell_all$sample_id[i] ] = tumor_cell_all$cor[i]
}

cell_info = unique(tumor_cell_all[,c("CCLE.name","tumor_type")])
rownames(cell_info) = cell_info$CCLE.name
row_annot = subset(cell_info, select= c("tumor_type"))
row_annot$tumor_type = as.character(row_annot$tumor_type)
row_annot$tumor_type[row_annot$tumor_type != cancer] = NA

pdf(paste(cancer, "/tumor_cell_line/", comparison_gene_set, "/tumor_cell_line.pdf", sep=""))
  
  reg1 = annHeatmap2(cell_tumor_matrix,               
                     ann=list(Row=list(data=row_annot[rownames(cell_tumor_matrix),])), 
                     scale="none",
                     cluster=list(Col=list(cuth=0.3)),
                     labels=list(Col = list(nrow=1,cex=0.01)),
                     legend=1)
  plot(reg1)
dev.off()

tumor_cell_mapping = read.csv("ccle_cellline_tcga_mapping_updated.csv")
row_annot = unique(subset(tumor_cell_mapping, select=c("CCLE.name", "Site.Primary", "Histology", "Hist.Subtype1", "tcga.tumor")))
rownames(row_annot) = row_annot$CCLE.name
row_annot$histology = "others"
row_annot$histology[row_annot$tcga.tumor == cancer] = cancer
row_annot$histology = as.factor(row_annot$histology)
row_annot = subset(row_annot, select=c(  "histology"))
#row_annot$tumor_type = as.factor(row_annot$tumor_type)
my.col = colorRampPalette(brewer.pal(8,"YlOrRd"))(100)
pheatmap(t(cell_tumor_matrix), annotation=row_annot, col= my.col, show_rownames=F, show_colnames=F, annotation_legend=T, file=paste(cancer, "/tumor_cell_line/", comparison_gene_set, "/tumor_cell_line_pheatmap.pdf", sep=""))
pheatmap(t(cell_tumor_matrix), annotation=row_annot, col= my.col, cellwidth=10, show_rownames=F, show_colnames=T, annotation_legend=T, file=paste(cancer, "/tumor_cell_line/", comparison_gene_set, "/tumor_cell_line_pheatmap_show_cell.pdf", sep=""))
pheatmap(t(cell_tumor_matrix), annotation=row_annot, col= my.col, cellwidth=10,cellheight=10, show_rownames=T, show_colnames=T, annotation_legend=T, file=paste(cancer, "/tumor_cell_line/", comparison_gene_set, "/tumor_cell_line_pheatmap_show_both.pdf", sep=""))


pdf(paste(cancer, "/tumor_cell_line/", comparison_gene_set, "/", cancer, "_all_celllines.pdf", sep=""))
  cell_line_type_order <- aggregate(cor ~ primary_site, tumor_cell_all, median)
  cell_line_type_ordered <- cell_line_type_order$primary_site[order(cell_line_type_order$cor, decreasing=T)]
  tumor_cell_all$primary_site = factor(tumor_cell_all$primary_site, levels = cell_line_type_ordered)
  p <- ggplot(tumor_cell_all, aes(primary_site, cor))
  print(p + geom_jitter(colour="black", size=0.5)  +  ylab("correlation") +
          stat_summary(geom = "crossbar", width=0.65, fatten=0, color="red", fun.data = function(x){ return(c(y=median(x), ymin=median(x), ymax=median(x))) }) + 
          xlab("") + theme(panel.background = element_rect(color = 'white'), axis.text.x = element_text(angle = 45, hjust = 1, size=15)) 
  )
dev.off()



#tumor_cell_all_subset = subset(tumor_cell_all, outlier < cutoff)
#limit to cancer specific cell lines
tumor_cell_all_subset = subset(tumor_cell_all, tumor_type == cancer & barcode %in% rownames(tumor_info[tumor_info$type == "tumor",]))

#test subset
cell_lines = unique(tumor_cell_all_subset$CCLE.name)
samples = unique(tumor_cell_all_subset$barcode)
cell_tumor_matrix = matrix(NA, nrow=length(cell_lines), ncol=length(samples), dimnames=list(cell_lines, samples))
for(i in 1:nrow(tumor_cell_all_subset)){
  cell_tumor_matrix[tumor_cell_all_subset$CCLE.name[i],tumor_cell_all_subset$barcode[i] ] = tumor_cell_all_subset$cor[i]
}


## 'Overwrite' default draw_colnames with your own version 
assignInNamespace(x="draw_colnames", value="draw_colnames_45",
                  ns=asNamespace("pheatmap"))

pheatmap(t(cell_tumor_matrix),  col= my.col, show_rownames=F, show_colnames=T, cex=0.6, cellwidth=8, file=paste(cancer, "/tumor_cell_line/", comparison_gene_set, "/tumor_cell_line_subset_pheatmap.pdf", sep=""))



pdf(paste(cancer, "/tumor_cell_line/", comparison_gene_set, "/", cancer, "_specific_celllines.pdf", sep=""))
  cell_line_order <- aggregate(cor ~ CCLE.name, tumor_cell_all_subset, median)
  cell_line_ordered <- cell_line_order$CCLE.name[order(cell_line_order$cor, decreasing=T)]
  tumor_cell_all_subset$CCLE.name = factor(tumor_cell_all_subset$CCLE.name, levels = cell_line_ordered)
  p <- ggplot(tumor_cell_all_subset, aes(CCLE.name, cor))
  print(p + geom_jitter(colour="black", size=1)  +  ylab("correlation to tumors") + theme_bw() + theme(axis.title.y = element_text(size = rel(1.5))) +
          stat_summary(geom = "crossbar", width=0.65, fatten=0, color="red", fun.data = function(x){ return(c(y=median(x), ymin=median(x), ymax=median(x))) }) + 
          xlab("") + theme(panel.background = element_rect(color = 'white'), axis.text.x = element_text(angle = 45, hjust = 1, size=10), axis.text.y = element_text(size=14)) 
  )
dev.off()

tumor_cell_all_subset$CCLE.name.truncated = sapply(as.character(tumor_cell_all_subset$CCLE.name), function(cell){
  unlist(strsplit(cell, "_"))[1]
})

pdf(paste(cancer, "/tumor_cell_line/", comparison_gene_set, "/", cancer, "_specific_celllines_violin.pdf", sep=""))
  cell_line_order <- aggregate(cor ~ CCLE.name.truncated, tumor_cell_all_subset, median)
  cell_line_ordered <- cell_line_order$CCLE.name[order(cell_line_order$cor, decreasing=T)]
  tumor_cell_all_subset$CCLE.name.truncated = factor(tumor_cell_all_subset$CCLE.name.truncated, levels = cell_line_ordered)
  p <- ggplot(tumor_cell_all_subset, aes(CCLE.name.truncated, cor))   
  print(p  + geom_violin( fill='grey',trim=F) +
          geom_jitter(colour="black",  size = 1, position=position_jitter(width=0.1, height=0)) +       
          ylab("correlation to tumors") + theme_bw() + theme(axis.title.y = element_text(size = rel(1.5))) +
          stat_summary(geom = "crossbar", width=0.65, fatten=0, color="red", fun.data = function(x){ return(c(y=median(x), ymin=median(x), ymax=median(x))) }) + 
          xlab("") + theme(panel.background = element_rect(color = 'white'), axis.text.x = element_text(angle = 45, hjust = 1, size=15), axis.text.y = element_text(size=16)) +
          theme(axis.title.y = element_text(size = rel(2), angle = 90))
  )
dev.off()

library(beanplot)
beanplot(cor ~ CCLE.name, data = tumor_cell_all_subset, 
         col = "lightgray", border = "grey", jitter=0.1, overallline = "median")
text(x = seq(0,1,length.out=length( tumor_cell_all_subset$CCLE.name ) ), c(-0.05),
     labels = c( unique(tumor_cell_all_subset$CCLE.name)), srt = 45, pos=2,offset=0.05, xpd = TRUE, cex=0.8)

cell_line_specific <- aggregate(cor ~ CCLE.name + Cell.line.primary.name, tumor_cell_all_subset, median)
cell_line_specific <- cell_line_specific[order(cell_line_specific$cor),]

all_cell_line <- aggregate(cor ~ CCLE.name + Cell.line.primary.name + tumor_type, tumor_cell_all, median)
cell_sig_cutoff <- qnorm(0.05, mean(all_cell_line$cor), sd(all_cell_line$cor), lower.tail = F)

cell_line_specific_sig <- subset(cell_line_specific, cor > cell_sig_cutoff)

write.csv(cell_line_specific_sig, paste(cancer, "/tumor_cell_line/", comparison_gene_set, "/", cancer, "_specific_cell_lines.csv", sep=""), row.names=F)
write.csv(cbind(all_cell_line, cutoff = cell_sig_cutoff), paste(cancer, "/tumor_cell_line/", comparison_gene_set, "/", cancer, "_all_cell_line.csv", sep=""), row.names=F)
