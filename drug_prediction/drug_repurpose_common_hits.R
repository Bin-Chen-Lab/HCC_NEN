#analyze predictions from CMap and LINCS
#drugbank should be updated (the version being used was downloaded in 2013)
##############################
######
library(pheatmap)
library(gplots)
library(ggplot2)
library("gplots")
library(RColorBrewer)
library("plyr")

#############################
#functions
############################
createGGPlotDist <- function(currPheno, measureVar, groupVars, main="GGPLot", ylab="pass ylab value", xlab="pass xlab value") {
  #modified from Purvesh's code
  
  temp = summarySE(currPheno, measurevar=measureVar, groupvars=groupVars)[, c(groupVars, measureVar, "se")]
  a = ggplot() +
    geom_violin(data=currPheno, aes(x=group,y=score), fill='grey',trim=F) +
    coord_flip() +
    geom_jitter(data=currPheno, aes(x=group,y=score, color=color), size = 2, position=position_jitter(width=jitterWidth, height=jitterHight)) +
    guides(color=FALSE) +
    xlab(xlab) +
    ylab(ylab) + 
    ggtitle(main) +  theme_bw() + 
    theme(text = element_text(size=24)) +
    #theme(legend.position=c(0.5,0.9), legend.title=element_blank(), legend.text=element_text(size=24)) +
    theme(axis.text.x = element_text(size=24, angle = 90, hjust = 1), axis.text.y = element_blank()) +
    geom_segment(data=temp,
                 aes(x=match(group,levels(group))-segmentWidth,
                     xend=match(group,levels(group))+segmentWidth,
                     y=score-se,yend=score-se),
                 col='black') +
    geom_segment(data=temp,
                 aes(x=match(group,levels(group))-segmentWidth,
                     xend=match(group,levels(group))+segmentWidth,
                     y=score+se,yend=score+se),
                 col='black') +
    geom_segment(data=temp,
                 aes(x=match(group,levels(group)),
                     xend=match(group,levels(group)),
                     y=score+se,yend=score-se),
                 col='black') +
    geom_point(data=temp, aes(x=group, y=score), color="black", size=1) 
    
  return (a)
}


createGGPlot <- function(currPheno, measureVar, groupVars, main="GGPLot", ylab="pass ylab value", xlab="pass xlab value", sig_cutoff=0) {
  temp = summarySE(currPheno, measurevar=measureVar, groupvars=groupVars)[, c(groupVars, measureVar, "se")]
  a = ggplot() +
    geom_violin(data=currPheno, aes(x=group,y=score), fill='white',trim=F) +
    coord_flip() +
    geom_jitter(data=currPheno, aes(x=group,y=score, color=cell_line), size = 2, position=position_jitter(width=jitterWidth, height=jitterHight)) +
    guides(shape=FALSE) +
    xlab(xlab) +
    ylab(ylab) + 
    ggtitle(main) +
    theme(text = element_text(size=24)) + theme_bw() +
    #theme(legend.position=c(0.5,0.9), legend.title=element_blank(), legend.text=element_text(size=24)) +
    theme(axis.text.x = element_text(size=24, angle = 90, hjust = 1), axis.text.y = element_text(size=24),  axis.title=element_text(size=24), legend.text = element_text(size=18), legend.title =element_text(size=18),  legend.position = c(0.8, 0.8)) +
    geom_segment(data=temp,
                 aes(x=match(group,levels(group))-segmentWidth,
                     xend=match(group,levels(group))+segmentWidth,
                     y=score-se,yend=score-se),
                 col='black') +
    geom_segment(data=temp,
                 aes(x=match(group,levels(group))-segmentWidth,
                     xend=match(group,levels(group))+segmentWidth,
                     y=score+se,yend=score+se),
                 col='black') +
    geom_segment(data=temp,
                 aes(x=match(group,levels(group)),
                     xend=match(group,levels(group)),
                     y=score+se,yend=score-se),
                 col='black') +
    geom_point(data=temp, aes(x=group, y=score), color="black", size=1) +
    geom_hline(yintercept = sig_cutoff, linetype = "longdash") + ylim(-1.5,1.5)
  return (a)
}


summarySE <- function(data=NULL, measurevar, groupvars=NULL, na.rm=FALSE,
                      conf.interval=.95, .drop=TRUE) {
  require(plyr)
  
  # New version of length which can handle NA's: if na.rm==T, don't count them
  length2 <- function (x, na.rm=FALSE) {
    if (na.rm) sum(!is.na(x))
    else       length(x)
  }
  
  # This does the summary. For each group's data frame, return a vector with
  # N, mean, and sd
  datac <- ddply(data, groupvars, .drop=.drop,
                 .fun = function(xx, col) {
                   c(N    = length2(xx[[col]], na.rm=na.rm),
                     mean = mean   (xx[[col]], na.rm=na.rm),
                     sd   = sd     (xx[[col]], na.rm=na.rm)
                   )
                 },
                 measurevar
  )
  
  # Rename the "mean" column    
  datac <- rename(datac, c("mean" = measurevar))
  
  datac$se <- datac$sd / sqrt(datac$N)  # Calculate standard error of the mean
  
  # Confidence interval multiplier for standard error
  # Calculate t-statistic for confidence interval: 
  # e.g., if conf.interval is .95, use .975 (above/below), and use df=N-1
  ciMult <- qt(conf.interval/2 + .5, datac$N-1)
  datac$ci <- datac$se * ciMult
  
  return(datac)
}

#########################
#MAIN
########################
cancer = "LIHC"

common_hits = read.csv(paste(cancer, "/drug/cmap_lincs_results_common.csv", sep=""))

hcc_drugs = read.csv("~/Documents/stanford/combination/drug/hcc_drugs.csv", stringsAsFactors = F)

common_hits[common_hits$x %in% tolower(hcc_drugs$drug),]
#############
###CMAP
load('raw/cmap/cmap_signatures.RData')
load(paste(cancer, "/drug/", "cmap_predictions.RData", sep=""))
valid_instances = read.csv("raw/cmap/cmap_valid_instances.csv", stringsAsFactors = F)
cmap_experiments = read.csv("raw/cmap/cmap_drug_experiments_new.csv", stringsAsFactors =  F)

drug_preds = results[[1]] #from cmap_predictions.RData
dz_signature = results[[2]]

cmap_experiments_valid = merge(cmap_experiments, valid_instances, by="id")
#keep valid profiles; drugs (as long as matched to Drugbank)
cmap_experiments_valid = cmap_experiments_valid[cmap_experiments_valid$valid == 1 & cmap_experiments_valid$DrugBank.ID != "NULL", ]
 
drug_instances_all = merge(drug_preds, cmap_experiments_valid, by.x="exp_id", by.y="id")
#write.csv(drug_instances_all, paste(cancer, "/drug/cmap_drug_predictions.csv", sep=""))

control = tail(drug_instances_all[order(drug_instances_all$cmap_score),], 1)$exp_id
#sig cmap cutoff
cutoff = max(drug_instances_all$cmap_score[drug_instances_all$q < 0.05 & drug_instances_all$cmap_score < 0])


drug_instances = subset(drug_instances_all, q < 0.05 & cmap_score < 0 & tolower(name) %in% common_hits$x )

drug_min_score = aggregate(cmap_score ~ name, drug_instances, min)

#
drug_instances = drug_instances[drug_instances$cmap_score %in% drug_min_score$cmap_score, ]
drug_instances = drug_instances[order(drug_instances$cmap_score), ]
drug_instances_id <- c(drug_instances$exp_id, control) + 1 #[c(1:20, (nrow(drug_instances)-5):nrow(drug_instances))] + 1 #the first column is the gene id
#get candidate drugs
drug_signatures <- cmap_signatures[,c(1, drug_instances_id)] #the first column is the gene id

drug_dz_signature <- merge(dz_signature[, c("GeneID", "value")], drug_signatures, by.x = "GeneID", by.y="V1")
drug_dz_signature <- drug_dz_signature[order(drug_dz_signature$value),]

#rerank 
drug_dz_signature[,2] = -drug_dz_signature[,2] # the higher rank indicate the more overexpressed, so we need to reverse the order
for (i in 2:ncol(drug_dz_signature)){
  drug_dz_signature[,i] = rank(drug_dz_signature[,i] )
}
drug_dz_signature <- drug_dz_signature[order(drug_dz_signature[,2]),] #order by disease expression

gene_ids <- drug_dz_signature[,1]
drug_dz_signature <- drug_dz_signature[, -1]

re_rank_instance = 0
if (re_rank_instance > 0){
  col_sorted = sort(cor(drug_dz_signature, method="spearman")["value",-1])    
  drug_dz_signature = drug_dz_signature[,c("value", names(col_sorted))]
}

drug_names <- sapply(2:ncol(drug_dz_signature), function(id){
  #need to minus 1 as in cmap_signatures, V1 is gene id.
  new_id = strtoi(paste(unlist(strsplit(as.character(colnames(drug_dz_signature)[id]),""))[-1], collapse="")) - 1 
  cmap_experiments_valid$name[cmap_experiments_valid$id == new_id]
})

in_vitro_list = c("niclosamide") #c("pyrvinium", "niclosamide", "mebendazole", "bithionol")
pdf(paste(cancer, "/drug/", "cmap_predictions_common.pdf", sep=""))
  colPal <- redgreen(100)
  par(mar=c(25, 4, 2, 0.5))
  #image(t(druggable_targets_pos), col=redblue(2))
  image(t(drug_dz_signature), col= colPal,   axes=F, srt=45)
  axis(1,  at=seq(0,1,length.out= ncol( drug_dz_signature ) ), labels= F)
  cols = ifelse(tolower(drug_names) %in% in_vitro_list, "red", "black")
  text(x = seq(0,1,length.out=ncol( drug_dz_signature ) ), c(-0.05),
     labels = c( "HCC",drug_names), srt = 45, pos=2,offset=0.05, xpd = TRUE, cex=1.4, col = c("red", cols))
dev.off()
write.csv(drug_dz_signature, paste(cancer, "/drug/cmap_drug_dz_signature.csv", sep=""))


################################################
##lincs
#####
###############################################
load('raw/lincs/lincs_signatures_cmpd_landmark.RData')
load(paste(cancer, "/drug/", "lincs_predictions.RData", sep=""))
landmark = read.csv("raw/lincs/lincs_landmark.csv", stringsAsFactors = F)
lincs_experiments = read.csv("raw/lincs/lincs_sig_info.csv", stringsAsFactors = F)

drug_preds = results[[1]]
dz_sig = results[[2]]

drug_preds$exp_id = as.numeric(as.character(drug_preds$exp_id))
drug_preds_sig = merge(drug_preds, lincs_experiments, by.x="exp_id", by.y="id")

lincs_predictions_all = subset(drug_preds_sig, cell_id %in% c("HEPG2", "HUH7") & DrugBank.ID != "NULL")
lincs_predictions_all = lincs_predictions_all[order(lincs_predictions_all$cmap_score),]

control = as.numeric(tail(lincs_predictions_all[order(lincs_predictions_all$cmap_score),], 1)$exp_id)


lincs_predictions = lincs_predictions_all
lincs_predictions = subset(lincs_predictions, q < 0.05 & cmap_score < 0 & tolower(pert_iname) %in% common_hits$x )
drug_min_score = aggregate(cmap_score ~ pert_iname, lincs_predictions, min)
lincs_predictions = lincs_predictions[lincs_predictions$cmap_score %in% drug_min_score$cmap_score, ]
lincs_predictions = lincs_predictions[order(lincs_predictions$cmap_score), ]
#visualize top 20 drugs
#get drug signatures
sig_ids = c(lincs_predictions$exp_id, control)  #[c(1:20, (nrow(lincs_predictions)-5):nrow(lincs_predictions))]
drug_names = lincs_predictions$pert_iname[1:20]
sigs = data.frame(gene_id =landmark$gene_id, lincs_signatures[, as.character(sig_ids)])

drug_dz_signature = merge(dz_sig[, c("GeneID", "value")], sigs, by.x = "GeneID", by.y="gene_id")
drug_dz_signature = drug_dz_signature[order(drug_dz_signature$value, decreasing = T),]
gene_ids = drug_dz_signature[,1]
drug_dz_signature = drug_dz_signature[, -c(1)]

#since the value is z score. the higher z score should have a higher rank...so we reverse the score
drug_dz_signature = -1 * drug_dz_signature
for (i in 1:ncol(drug_dz_signature)){
  drug_dz_signature[,i] = rank(drug_dz_signature[,i])
}

re_rank_instance = 0
if (re_rank_instance > 0){
  col_sorted = sort(cor(drug_dz_signature, method="spearman")["value",-1])    
  drug_dz_signature = drug_dz_signature[,c("value", names(col_sorted))]
}

drug_names <- sapply(2:ncol(drug_dz_signature), function(id){
  new_id = strtoi(paste(unlist(strsplit(as.character(colnames(drug_dz_signature)[id]),""))[-1], collapse="")) 
  lincs_experiments$pert_iname[lincs_experiments$id == new_id]
})


gene_annot <- subset(dz_sig, select = c("Symbol", "GeneID"))

in_vitro_list = c("niclosamide") #c("pyrvinium", "niclosamide", "mebendazole", "bithionol")
###visualze top 20 genes
pdf(paste(cancer, "/drug/lincs_predictions_common.pdf", sep=""))
  layout(matrix(1))
  par(mar=c(25, 4, 1, 0.5))
  colPal <- redgreen(100)
  image(t(drug_dz_signature), col= colPal,   axes=F, srt=45)
  cols = ifelse(tolower(drug_names) %in% in_vitro_list, "red", "black")
  axis(1,  at=seq(0,1,length.out= ncol( drug_dz_signature ) ), labels= F)
  text(x = seq(0,1,length.out=ncol( drug_dz_signature ) ), c(-0.05),
       labels = c( "HCC",drug_names), srt = 45, pos=2,offset=0.05, xpd = TRUE, cex=1.4, col = c("red",   cols = ifelse(tolower(drug_names) %in% in_vitro_list, "red", "black")))
dev.off()

