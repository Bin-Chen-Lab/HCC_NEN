#prepare signatures for drug predictions

cancer <- "LIHC"
load(paste(cancer, "/", cancer, '_sig_deseq_final.RData', sep=''))
res$symbol= sapply(res$id, function(id){
  as.character(unlist(strsplit(id, "\\|"))[1])
})
res$GeneID <- sapply(res$id, function(id){
  as.character(unlist(strsplit(id, "\\|"))[2])
})

if (!file.exists(paste(cancer, "/drug", sep=""))){
  dir.create(paste(cancer, "/drug", sep=""))
}

#for cmap, more stringent
dz_signature <- subset(res, !is.na(padj) & !is.na(id) & id !='?' & padj < 1E-12 & abs(log2FoldChange) > 2 & abs(log2FoldChange) != Inf)

dz_signature <- subset(dz_signature, select=c("GeneID","symbol", "log2FoldChange", "padj"))
names(dz_signature) <- c("GeneID", "Symbol", "value", "padj")
dz_signature <- subset(dz_signature, !is.na(value))
dz_signature <- dz_signature[order(dz_signature$value),]

#emperical experience. Signatures around ~300 are suitable for cmap
if (nrow(dz_signature) > 300){
  dz_signature <- rbind(head(dz_signature, 150), tail(dz_signature, 150))
}

dz_signature$up_down <- "up"
dz_signature$up_down[dz_signature$value < 0] <- "down"
write.table(dz_signature, paste(cancer, "/drug/","dz_signature_cmap.txt", sep=""), quote=F, col.names=T, row.names=F, sep="\t" )

############
#for lincs, 
dz_signature <- subset(res, !is.na(padj) & !is.na(id) & id !='?' & padj < 1E-03 & abs(log2FoldChange) > 2 & abs(log2FoldChange) != Inf)

dz_signature <- subset(dz_signature, select=c("GeneID","symbol", "log2FoldChange"))
names(dz_signature) <- c("GeneID", "Symbol", "value")
dz_signature <- subset(dz_signature, !is.na(value))
dz_signature <- dz_signature[order(dz_signature$value),]
dz_signature$up_down <- "up"
dz_signature$up_down[dz_signature$value < 0] <- "down"

landmark <- read.csv("raw/lincs//lincs_landmark.csv")
dz_signature <- subset(dz_signature, GeneID %in% landmark$gene_id)
write.table(dz_signature, paste(cancer, "/drug/","dz_signature_lincs.txt", sep=""), quote=F, col.names=T, row.names=F, sep="\t" )
