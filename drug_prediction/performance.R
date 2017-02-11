#setwd("~/Documents/stanford/hcc/doc/niclosamide_manuscript/doc/hepatology/supp/")
#enrichment of known drugs

library("ROCR")
library("GSVA")
library("colorRamps")
hcc_drugs = read.csv("~/Documents/stanford/combination/drug/hcc_drugs.csv", stringsAsFactors = F)

database = "cmap"

if (database == "cmap"){
#CMap results
load("~/Documents/stanford/hcc/release/data/LIHC/drug/cmap_predictions.RData")
drug_preds = results[[1]]
dz_sig = results[[2]]
cmap_experiments = read.csv("~/Documents/stanford/hcc/release/data/raw/cmap/cmap_drug_experiments.csv")
valid_instances = read.csv("~/Documents/stanford/hcc/release/data/raw/cmap/cmap_valid_instances.csv")
#cmap_experiments = cmap_experiments[valid_instances$valid == 1, ]
cmap_prediction = merge(drug_preds, cmap_experiments, by.x="exp_id", by.y="id")
cmap_prediction = subset(cmap_prediction,  DrugBank.ID != "NULL") #&     & DrugBank.ID != "NULL"

}else{
#LINCS
load("~/Documents/stanford/hcc/release/data/LIHC/drug/lincs_predictions.RData")
landmark = read.csv("~/Documents/stanford/hcc/release/data/raw/lincs/lincs_landmark.csv", stringsAsFactors = F)
lincs_experiments = read.csv("~/Documents/stanford/hcc/release/data/raw/lincs/lincs_sig_info.csv", stringsAsFactors = F)
drug_preds = results[[1]]
dz_sig = results[[2]]
cmap_prediction = merge(drug_preds, lincs_experiments, by.x="exp_id", by.y="id")
cmap_prediction = subset(cmap_prediction, cell_id %in% c("HEPG2", "HUH7") & DrugBank.ID != "NULL") #&    &   & DrugBank.ID != "NULL"
cmap_prediction$name = cmap_prediction$pert_iname
}

cmap_prediction = aggregate(cmap_score ~ name, cmap_prediction, min)
cmap_prediction$drug = 0
cmap_prediction$drug[tolower(cmap_prediction$name) %in% tolower(hcc_drugs$drug)] = 1

pred.roc=prediction(cmap_prediction$cmap_score,cmap_prediction$drug)
performance(pred.roc,"auc")


cmap_prediction = cmap_prediction[order(cmap_prediction$cmap_score), ]
cmap_prediction[cmap_prediction$drug == 1, ]

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

par(mar=c(12, 1, 1, 1))

#adjust cmap score for better visualization

cmap_prediction$score = sapply(cmap_prediction$cmap_score, function(x){
  if (x < 0){
    -20^(abs(x))
  }else{
    20^(abs(x))
  }
})
pdf( paste( "fig/enrichment_", database, "_hcc_drugs.pdf", sep=""))
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

pdf( paste( "fig/enrichment_", database, ".pdf", sep=""))
par(mar=c(12, 4, 12, 4)) #bottom, left, top, and right
z = (matrix(cmap_prediction$score))
n = 100
image(x=1:nrow(z), y = 1, z, col = green2red(n), xlab = "", ylab="", axes=FALSE )

#text(-50, 1, "aa")
dev.off()
