#select disease signatures (p adj and fc) that lead to the best separation between tumors and non-tumors.

library(pheatmap)
library(scatterplot3d) 
library(pROC)
library(preprocessCore) # include quantile normalization
library(ggplot2)

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
  if (is.null(ct$Logged)) {ct$Logged = NA}
  
  return(ct)
}

####################
#MAIN
###################

###
#disease signatures
cancer <- "LIHC"
dz_padj_cutoffs <-  10^-seq(2, 20)
dz_fc_cutoffs <- c(1, 1.5, 2, 2.5, 3)

load(paste(cancer, "/", cancer, '_sig_deseq_final.RData', sep=''))
res$symbol= sapply(res$id, function(id){
  as.character(unlist(strsplit(id, "\\|"))[1])
})
res$GeneID <- sapply(res$id, function(id){
  as.character(unlist(strsplit(id, "\\|"))[2])
})


ct<-readCodingTable("raw/geo/meta_input.txt")
ct <- subset(ct, valid == 1)

OriginCodes <- unique(ct$OriginCode)

performances <- data.frame()
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
  
  for (dz_padj_cutoff in dz_padj_cutoffs){
    for (dz_fc_cutoff in dz_fc_cutoffs ){
      dz_signature <- subset(res, !is.na(padj) & !is.na(id) & id !='?' & padj < dz_padj_cutoff & abs(log2FoldChange) > dz_fc_cutoff & abs(log2FoldChange) != Inf )
      
      expressions_subset <- expressions[expressions$GeneID %in% dz_signature$GeneID, colnames(expressions) %in% c("GeneID", ct_subset$GSM)]
      
      expressions_by_gene <- aggregate(. ~ GeneID, data= expressions_subset, mean)
      geneids <- expressions_by_gene$GeneID
      expressions_by_gene <- expressions_by_gene[, -1]
      genenames <- merge(geneids, dz_signature, by.x=1, by.y="GeneID", sort=F)$Symbol
      rownames(expressions_by_gene) <- genenames
      
      annotation <- subset(ct_subset, select= c("GSM", "ClassCode"))
      rownames(annotation) <- annotation$GSM
      annotation$ClassCode <- as.factor(annotation$ClassCode)
      annotation <- subset(annotation, select=c("ClassCode"))
      
      pca <- prcomp(t(expressions_by_gene))
      sample_class <- annotation[rownames(pca$x),1]
      roc_result <- roc(as.character(sample_class), scale(as.numeric(pca$x[,1])) )
      performances <- rbind(performances, data.frame(code,dz_padj_cutoff, dz_fc_cutoff, roc = roc_result$auc, size = length(dz_signature$GeneID)))
    }
  }
}

if (!file.exists(paste(cancer, "/dz_sig", sep=""))){
  dir.create(paste(cancer, "/dz_sig", sep=""))
}

write.csv(performances,paste(cancer, "/dz_sig/", "performances.csv", sep=""))

#performances = read.csv(paste(cancer, "/dz_sig/", "performances.csv", sep=""))

#this shows fc = 2 performs the best
pdf(paste(cancer, "/dz_sig/", "fc_auc.pdf", sep=""))
  performances$dataset <- performances$code
  performances$cutoff_code <- paste(performances$dataset, performances$dz_padj_cutoff)
  h <- ggplot(performances, aes(dz_fc_cutoff, roc))
  h + geom_line(aes(group =  cutoff_code, color = dataset)) + xlab("log2 fold change") + ylab("AUC")
dev.off()

#this shows padj doe not have any effect
pdf(paste(cancer, "/dz_sig/","p_auc.pdf", sep=""))
  performances$cutoff_code <- paste(performances$dataset, performances$dz_fc_cutoff)
  h <- ggplot(performances, aes(log(dz_padj_cutoff, 10), roc))
  h + geom_line(aes(group =  cutoff_code, color = dataset)) + xlab("log10 p adjust") + ylab("AUC")
dev.off()

log(dz_padj_cutoff, 10)
a <- subset(performances, dz_padj_cutoff == 1e-10)
plot((a$dz_fc_cutoff), a$roc)

a <- subset(performances, dz_fc_cutoff == 2)
plot(log(a$dz_padj_cutoff, 10), a$roc)
a_order <- aggregate(roc ~ dz_padj_cutoff, a, median)
a_order <- a_order[order(a_order$roc), ]
tail(a_order)

pdf(paste(cancer, "/dz_sig/", "p_auc_median.pdf", sep=""))
  h <- ggplot(a_order, aes(log(a_order$dz_padj_cutoff, 10), roc))
  h + geom_line() + geom_point() + xlab("log10 p adjust") + ylab("AUC")
dev.off()
