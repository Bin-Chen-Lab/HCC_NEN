#we find that many intances of one drug behavior quite differently, we should ignore the instances inconsistent with others
find_validated_sample_by_cor <- function(results,min.cor){
  if (ncol(results)<2) return(NA)
  sim <- cor(results,method="spearman")
  diag(sim) <- 0
  sig_samples <- sapply(1:nrow(sim),function(id){
    if (max(sim[id,]) >= min.cor) row.names(sim)[id]
  })
  sig_samples <- unlist(sig_samples)  
}

get.significant.cor.cutoff <- function(pairs_cor,pvalue){
  #given a correlation matrix, find the quantile greater than p value
  meanCorr=apply(pairs_cor,2,mean)  
  sdCorr = apply(pairs_cor,2,sd)
  
  return(qnorm(pvalue,mean(meanCorr),mean(sdCorr),lower.tail=F))
}

get.drug.instance <- function(id,min.cor){
  return(id)
}

load('~/Documents/stanford/lincs/data/cmap_signatures_updated.RData')
ls()
cmap_signatures[1:2,1:2]
gene_ids <- cmap_signatures[,1]
cmap_signatures <- cmap_signatures[,-1]
all_instances_cor <- cor(cmap_signatures,method="spearman")


q <- get.significant.cor.cutoff(all_instances_cor,0.05)

#get instance info
library("RMySQL")
mysql_drvr <-dbDriver("MySQL")
#con <- dbConnect(mysql_drvr,user="binchen1",password="binchenstanford",host="buttelab-db1.stanford.edu",dbname="user_binchen1")

rs <- dbSendQuery(con, "select distinct cmap_id as instance_id, pert_desc as drug from v_lincs_cmap_id_info")
lincs <- fetch(rs, n = -1)

rs <- dbSendQuery(con, "select distinct id as instance_id, name as drug from proj_repositioning.cmap_experiments")
cmap <- fetch(rs, n = -1)

cmap_lincs <- rbind(cmap,lincs)

drugs = unique(cmap_lincs$drug)
valid_instances = NULL
for (drug in drugs){
  print(drug)
  instances = cmap_lincs$instance_id[cmap_lincs$drug==drug]
  instance_sigs = subset(cmap_signatures,select=instances)
  valid_instance = find_validated_sample_by_cor(instance_sigs,q)
  #at least two instance per drug
  if (length(valid_instance)>1) {
    valid_instances = c(valid_instances,valid_instance)    
  }
}

cmap_lincs_valid_instances = data.frame(id = seq(1,ncol(cmap_signatures)),name =colnames(cmap_signatures),valid=(colnames(cmap_signatures) %in% valid_instances))
#dbWriteTable(con,"proj_repositioning.cmap_lincs_valid_instances",cmap_lincs_valid_instances)
