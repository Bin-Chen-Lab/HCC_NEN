#setwd("~/Documents/stanford/hcc/doc/niclosamide_manuscript/doc/hepatology/supp/")
#enrichment of known drugs

library("ROCR")
library("GSVA")
library("colorRamps")
hcc_drugs = read.csv("~/Documents/stanford/combination/drug/hcc_drugs.csv", stringsAsFactors = F)

database = "lincs"
landmark_used = 0
cancer = "LIHC"
if (database == "cmap"){
#predictions
  load('raw/cmap/cmap_signatures.RData')
  load(paste(cancer, "/drug/", "cmap_predictions.RData", sep=""))
  valid_instances = read.csv("raw/cmap/cmap_valid_instances.csv", stringsAsFactors = F)
  cmap_experiments = read.csv("raw/cmap/cmap_drug_experiments_new.csv", stringsAsFactors =  F)
  
  drug_preds = results[[1]] #from cmap_predictions.RData
  
  dz_signature = results[[2]]
  dz_signature = dz_signature[, c("GeneID", "value")]
  
  dz_cmap = merge(dz_signature, cmap_signatures, by.x="GeneID", by.y = "V1" )
  
  drug_preds_other = cor(-1 * dz_cmap$value, dz_cmap[, 3:ncol(dz_cmap)], method="spearman") #drug signature is absolute value
  
  drug_preds_other = data.frame(exp_id = 1:6100, cmap_score = as.numeric(drug_preds_other))


  cmap_experiments = read.csv("~/Documents/stanford/hcc/release/data/raw/cmap/cmap_drug_experiments.csv")
  valid_instances = read.csv("~/Documents/stanford/hcc/release/data/raw/cmap/cmap_valid_instances.csv")
  #cmap_experiments = cmap_experiments[valid_instances$valid == 1, ]
  cmap_prediction = merge(drug_preds_other, cmap_experiments, by.x="exp_id", by.y="id")
  cmap_prediction = subset(cmap_prediction,  DrugBank.ID != "NULL") #&     & DrugBank.ID != "NULL"
  #cmap_prediction = cmap_prediction[cmap_prediction$valid ==1, ]

}else{
#LINCS
  if (landmark_used == 1){
    load('raw/lincs/lincs_signatures_cmpd_landmark.RData')
    load("~/Documents/stanford/hcc/release/data/LIHC/drug/lincs_predictions.RData")
  }else{
    load('raw/lincs/lincs_signatures_cmpd_hcc_all.RData')
    load("~/Documents/stanford/hcc/release/data/LIHC/drug/cmap_predictions.RData") #need it to retrieve disease signature
  }
  
  landmark = read.csv("~/Documents/stanford/hcc/release/data/raw/lincs/lincs_landmark.csv", stringsAsFactors = F)
  lincs_experiments = read.csv("~/Documents/stanford/hcc/release/data/raw/lincs/lincs_sig_info.csv", stringsAsFactors = F)

  dz_sig = results[[2]]
  
  dz_signature = results[[2]]
  dz_signature = dz_signature[, c("GeneID", "value")]
  
  dz_cmap = merge(dz_signature, data.frame(GeneID = rownames(lincs_signatures), lincs_signatures), by ="GeneID" )
  colnames(dz_cmap) = c("GeneID", "value", colnames(lincs_signatures))
  
  if (landmark_used == 1){
    drug_preds_other = cor(1 * dz_cmap$value, dz_cmap[, 3:ncol(dz_cmap)], method="spearman") #drug signature is absolute value
  }else{
    drug_preds_other = cor(-1 * dz_cmap$value, dz_cmap[, 3:ncol(dz_cmap)], method="spearman") #drug signature is ranked
  }
  drug_preds_other = data.frame(exp_id =  colnames(drug_preds_other), cmap_score = as.numeric(drug_preds_other))
  
  cmap_prediction = merge(drug_preds_other, lincs_experiments, by.x="exp_id", by.y="id")
  cmap_prediction = subset(cmap_prediction, cell_id %in% c("HEPG2", "HUH7") & DrugBank.ID != "NULL") #&    &   & DrugBank.ID != "NULL"
  cmap_prediction$name = cmap_prediction$pert_iname
}

cmap_prediction = aggregate(cmap_score ~ name, cmap_prediction, mean)
cmap_prediction$drug = 0
cmap_prediction$drug[tolower(cmap_prediction$name) %in% tolower(hcc_drugs$drug)] = 1

pred.roc=prediction(cmap_prediction$cmap_score,cmap_prediction$drug)

cmap_prediction[cmap_prediction$drug == 1, ]

cmap_prediction = cmap_prediction[order(cmap_prediction$cmap_score), ]

geneSets = list(hcc = cmap_prediction[cmap_prediction$drug == 1, "name"])

times = 10000
ranks = matrix(NA, nrow = nrow(cmap_prediction), ncol = times)
ranks[,1] = rank(cmap_prediction$cmap_score)
rownames(ranks) = cmap_prediction$name

for (i in 2:times){
  ranks[,i] = rank(sample(1:length(cmap_prediction$cmap_score), length(cmap_prediction$cmap_score)))
}

gsea_results = gsva(ranks, geneSets, method = "ssgsea")

sum(gsea_results[1, -1] < gsea_results[1,1])/length(gsea_results[1, -1])

'
par(mar=c(12, 1, 1, 1))

#adjust cmap score for better visualization

cmap_prediction$score = sapply(cmap_prediction$cmap_score, function(x){
  if (x < 0){
    -20^(abs(x))
  }else{
    20^(abs(x))
  }
})
pdf( paste( "fig/enrichment_", database, "_hcc_drugs_spearman.pdf", sep=""))
  par(mar=c(12, 4, 12, 4)) #bottom, left, top, and right
  z = (matrix((cmap_prediction$score)))
  n = 200
  image(x=1:nrow(z), y = 1, z, col = green2red(n), xlab = "", ylab="", axes=FALSE )
  drugs_pos = which(cmap_prediction$drug == 1)
  for (drug_pos in drugs_pos){
    abline(v = drug_pos)
  }
#text(-50, 1, "aa")
dev.off()

pdf( paste( "fig/enrichment_", database, "_spearman.pdf", sep=""))
  par(mar=c(12, 4, 12, 4)) #bottom, left, top, and right
  z = (matrix(cmap_prediction$score))
  n = 100
  image(x=1:nrow(z), y = 1, z, col = green2red(n), xlab = "", ylab="", axes=FALSE )

#text(-50, 1, "aa")
dev.off()
'