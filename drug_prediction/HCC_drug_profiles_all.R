#retrieve all HCC gene expression
get.instance.sig <- function(id, con,  landmark=F){
  
  sig_file = paste("~/Documents/stanford/lincs/data/lincs/", "id", ".txt", sep="")
  sig_value = NULL
  
  if (file.exists(sig_file)){
    sig_value= scan(sig_file )
    
    if (landmark){
      sig_value = sig_value[1:978]
    }
    
  }else{
    query <- paste("select * from proj_lincs.sig_values where id = ", id, sep="")
    rs <- dbSendQuery(con, query)
    result <- fetch(rs, n = -1)[1,]
    value <- as.double(unlist(strsplit(result[1,2], ",")))
    if (landmark == T){
      instance.sig <- data.frame(sig_id =  id, probe_id = seq(1, 978), value[1:978])      
    }else{
      instance.sig <- data.frame(sig_id = id, probe_id = seq(1, length(value)), value = value)      
    }
    dbClearResult(rs)  
    sig_value = instance.sig$value
  }
  return (sig_value)
}

get.gene.list <- function(con){
  lincs_gene <- dbReadTable(con, "probe_id_info")
  
  return(lincs_gene$gene_id)
}

lincs_experiments = read.csv("raw/lincs/lincs_sig_info.csv", stringsAsFactors = F)
lincs_experiments = lincs_experiments[lincs_experiments$cell_id %in% c("HEPG2", "HUH7"),]

exp_ids = lincs_experiments$id

cmap_exp_signature = data.frame(GeneID = get.gene.list(con))

count = 0 
for (exp_id in exp_ids){
  count = count + 1
  print(count)
  sig <-  get.instance.sig(exp_id, con)
  cmap_exp_signature <- cbind(cmap_exp_signature,  x=rank(-1 * sig, ties.method="random"))
}
colnames(cmap_exp_signature) = c("GeneID",exp_ids )

lincs_signatures = cmap_exp_signature[,-1]
rownames(lincs_signatures) = lincs_signatures[,1]
save(lincs_signatures, file = 'raw/lincs/lincs_signatures_cmpd_hcc_all.RData')

