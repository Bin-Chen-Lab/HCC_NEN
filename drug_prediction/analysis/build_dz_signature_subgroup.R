setwd("/Users/binchen1/Documents/stanford/hcc/data")
library("RMySQL")
mysql_drvr <-dbDriver("MySQL")
#con <- dbConnect(mysql_drvr,group="client",host="buttelab-db1.stanford.edu",dbname="proj_lincs")
get.landmark.info <- function(con){
  rs <- dbSendQuery(con, "select * from probe_id_info where pool_id like '%epsilon%'")    
  probe_info <- fetch(rs, n = -1)
  dbClearResult(rs)
  return(probe_info)
}
lincs_landmark = get.landmark.info(con)


cancer = "LIHC"
subtype = 'G1'
#load(paste( cancer, "/drug_subgroup/", subtype, "/Alcohol consumption_sig_deseq_matched.RData", sep=''))
#'Hepatitis B', 'Hepatitis C', 'Alcohol consumption'
#load(paste("sigs", '/', subtype , '_sig_deseq_final.RData', sep='')) #tcga_deseq LIHC_sig_deseq_QC
load(paste(cancer, '/', subtype , '_sig_deseq_final_all_non_tumors.RData', sep='')) #tcga_deseq LIHC_sig_deseq_QC


res$GeneID = sapply(res$id, function(id){
  unlist(strsplit(id, "\\|"))[2]
})
res$symbol = sapply(res$id, function(id){
  unlist(strsplit(id, "\\|"))[1]
})



#for cmap, more stringent
dz_signature = subset(res, !is.na(padj) & !is.na(id) & id !='?' & padj < 1E-2 & abs(log2FoldChange) > 1.5 & abs(log2FoldChange) != Inf)
#dz_signature = subset(dz_signature, baseMean > quantile(res$baseMean)[2])
dim(dz_signature)

dz_signature = subset(dz_signature, select=c("GeneID","symbol", "log2FoldChange"))
names(dz_signature) = c("GeneID", "Symbol", "value")
dz_signature = subset(dz_signature, !is.na(value))
dz_signature = dz_signature[order(dz_signature$value),]
if (nrow(dz_signature) > 300){
  dz_signature = rbind(head(dz_signature, 150), tail(dz_signature, 150))
}

dz_signature$up_down = "up"
dz_signature$up_down[dz_signature$value < 0] = "down"
write.table(dz_signature, paste(cancer, "/drug/",cancer, subtype, "_dz_signature_cmap_final_all_tumor_15_1E-2.txt", sep=""), quote=F, col.names=T, row.names=F, sep="\t" )

#for lincs, more loose, 
dz_signature = subset(res, !is.na(padj) & !is.na(id) & id !='?' & padj < 1E-2 & abs(log2FoldChange) > 1.5 & abs(log2FoldChange) != Inf)
#dz_signature = subset(dz_signature, baseMean > quantile(res$baseMean)[2])

dz_signature = subset(dz_signature, select=c("GeneID","symbol", "log2FoldChange"))
names(dz_signature) = c("GeneID", "Symbol", "value")
dz_signature = subset(dz_signature, !is.na(value))
dz_signature = dz_signature[order(dz_signature$value),]
dz_signature$up_down = "up"
dz_signature$up_down[dz_signature$value >  0] = "down"
dz_signature = subset(dz_signature, GeneID %in% lincs_landmark$gene_id)

write.table(dz_signature, paste(cancer, "/drug/",cancer, subtype, "_dz_signature_lincs_final_all_tumor_15_1E-2.txt", sep=""), quote=F, col.names=T, row.names=F, sep="\t" )
