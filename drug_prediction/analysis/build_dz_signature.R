setwd("/Users/binchen1/Documents/stanford/hcc/data")

cancer = "LIHC"
load(paste('sigs/', cancer, '_sig_deseq_final.RData', sep=''))
res$symbol= sapply(res$id, function(id){
  as.character(unlist(strsplit(id, "\\|"))[1])
})
res$GeneID = sapply(res$id, function(id){
  as.character(unlist(strsplit(id, "\\|"))[2])
})


#for cmap, more stringent
dz_signature = subset(res, !is.na(padj) & !is.na(id) & id !='?' & padj < 1E-12 & abs(log2FoldChange) > 2 & abs(log2FoldChange) != Inf)
#dz_signature = subset(dz_signature, baseMean > quantile(res$baseMean)[2])

dz_signature = subset(dz_signature, select=c("GeneID","symbol", "log2FoldChange", "padj"))
names(dz_signature) = c("GeneID", "Symbol", "value", "padj")
dz_signature = subset(dz_signature, !is.na(value))
dz_signature = dz_signature[order(dz_signature$value),]
if (nrow(dz_signature) > 300){
  dz_signature = rbind(head(dz_signature, 150), tail(dz_signature, 150))
}

dz_signature$up_down = "up"
dz_signature$up_down[dz_signature$value < 0] = "down"
write.table(dz_signature, paste(cancer, "/drug/",cancer, "_dz_signature_cmap_final_v1.txt", sep=""), quote=F, col.names=T, row.names=F, sep="\t" )

#for lincs, more loose, 
#cutoffs leading to the best median AUC of 6 datasets are: FC: 1, q = 1E-05, auc: 0.989
dz_signature = subset(res, !is.na(padj) & !is.na(id) & id !='?' & padj < 1E-05 & abs(log2FoldChange) > 1 & abs(log2FoldChange) != Inf)
#dz_signature = subset(dz_signature, baseMean > quantile(res$baseMean)[2])

dz_signature = subset(dz_signature, select=c("GeneID","symbol", "log2FoldChange"))
names(dz_signature) = c("GeneID", "Symbol", "value")
dz_signature = subset(dz_signature, !is.na(value))
dz_signature = dz_signature[order(dz_signature$value),]
dz_signature$up_down = "up"
dz_signature$up_down[dz_signature$value < 0] = "down"

landmark = read.csv("LIHC/moa/lincs_landmark.csv")
dz_signature = subset(dz_signature, GeneID %in% landmark$gene_id)
write.table(dz_signature, paste(cancer, "/drug/",cancer, "_dz_signature_lincs_1_1E-05.txt", sep=""), quote=F, col.names=T, row.names=F, sep="\t" )
