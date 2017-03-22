#create ROC curves for external datasets using disease signatures 

library(pheatmap)
library(scatterplot3d) 
library(pROC)

###################
#functions
###################
readCodingTable <- function(infile, DEBUG=0) {
  # Coding table needs to have a ClassCode (1 or 2), an Origin Code (some unique term for its identifier), GSM (particular sample ID), and GPL (the GEO platform ID)
  ct <- read.csv(infile, sep="\t", as.is=T)
  if (is.null(ct$ClassCode) | is.null(ct$OriginCode) | is.null(ct$GSM) | is.null(ct$GPL) )  {
    print("ERROR!!  Problem with coding table, needs to have columns for ClassCode (1/2 for treatment/control or a zero, negative or NA for a two color), OriginCode, GSM, and GPL.  This information is all essential")
    return(NA)		
  }
  if (is.null(ct$Logged)) {ct$Logged <- NA}
  
  return(ct)
}

####################
#MAIN
###################


###
#disease signatures
cancer <- "LIHC"
dz_padj_cutoff <- 1E-12
dz_fc_cutoff <- 2
load(paste(cancer, "/", cancer, '_sig_deseq_final.RData', sep=''))
res$symbol= sapply(res$id, function(id){
  as.character(unlist(strsplit(id, "\\|"))[1])
})
res$GeneID <- sapply(res$id, function(id){
  as.character(unlist(strsplit(id, "\\|"))[2])
})
dz_signature <- subset(res, !is.na(padj) & !is.na(id) & id !='?' & padj < dz_padj_cutoff & abs(log2FoldChange) > dz_fc_cutoff & abs(log2FoldChange) != Inf )


ct<-readCodingTable("meta_input.txt")
ct <- subset(ct, valid == 1)

OriginCodes <- unique(ct$OriginCode)

rocs <- list()
  
for (code in OriginCodes){
  ct_subset <- subset(ct, OriginCode == code)
  #find GSE
  GSE <- ct_subset$GSE[1]
  #find GPL
  GPL <- ct_subset$GPL[1]
  #map
  
  #Pull out raw expression
  rawValues <- read.csv(paste("raw/geo/rawExpression/", GPL, ".RawExpressionValues.table", sep=""), sep="\t" )
  raw_matrix1 <- normalize.quantiles(as.matrix(rawValues ))
  row.names(raw_matrix1) <- row.names(rawValues)
  colnames(raw_matrix1) <- colnames(rawValues)
  rawValues <- raw_matrix1
  
  #pull out gene mapping
  geneMappings <- read.csv(paste("raw/geo/geneMappings/", GPL, ".Probe2EntrezMap.table", sep=""), sep="\t" )
  geneMappings <- subset(geneMappings, select=c("probe", "GeneID"))
  #build matrix
  
  rawValues <- data.frame(probe_id=rownames(rawValues), rawValues)
  
  expressions <- merge(rawValues, geneMappings, by.x="probe_id", by.y="probe")
  
  expressions <- expressions[expressions$GeneID %in% dz_signature$GeneID, colnames(expressions) %in% c("GeneID", ct_subset$GSM)]
  
  expressions_by_gene <- aggregate(. ~ GeneID, data= expressions, mean)
  geneids <- expressions_by_gene$GeneID
  expressions_by_gene <- expressions_by_gene[, -1]
  genenames <- merge(geneids, dz_signature, by.x=1, by.y="GeneID", sort=F)$Symbol
  rownames(expressions_by_gene) <- genenames
  
  annotation <- subset(ct_subset, select= c("GSM", "ClassCode"))
  rownames(annotation) <- annotation$GSM
  annotation$ClassCode <- as.factor(annotation$ClassCode)
  annotation <- subset(annotation, select=c("ClassCode"))
  
  my.cols <- greenred(100) # brewer.pal(9, "Blues")
 # pheatmap(scale(expressions_by_gene), col = my.cols, annotation = annotation,  cellheight= 12, show_colnames=F, legend=F, 
 #          show_rownames=F , filename=paste( cancer,"/dz_sig/heatmap_", code, ".pdf", sep="")) #
  
  #if (nrow(expressions_by_gene) < ncol(expressions_by_gene)){
    pca <- prcomp(t(expressions_by_gene))
    sample_class <- annotation[rownames(pca$x),1]
    
    pdf(paste( cancer,"/dz_sig/pca_", code, "_final.pdf", sep=""))
    scatterplot3d(pca$x[,1:3], pch=20, color=sample_class) 
    dev.off()
  #}
  
    rocs[[code]] <- data.frame(pred = pca$x[,1], label = sample_class)
}

save(rocs, file = paste(cancer, "/dz_sig/rocs.RData", sep=""))

#load( paste(cancer, "/dz_sig/rocs.RData", sep=""))

pdf(paste(cancer, "/dz_sig/rocs.pdf", sep=""))
  plot(roc(label ~ pred, rocs[[1]],col=1))
  
  for (i in 2:length(rocs)){
    plot(roc(label ~ pred, rocs[[i]]), add=T, col=i)
  }
  legend("bottomright", c(names(rocs)), lty=1,  col=c(1:length(rocs)), cex=1.4)
dev.off()


aucs <- sapply(1:length(rocs), function(i){
  roc <- roc(label ~ pred, rocs[[i]],col=1)
  roc$auc
})


mean(aucs)
#0.96
#0.98
#0.995

