#analyze predictions from CMap and LINCS

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

#############
###CMAP
load('raw/cmap/cmap_signatures.RData')
load(paste(cancer, "/drug/", "cmap_predictions.RData", sep=""))
valid_instances = read.csv("raw/cmap/cmap_valid_instances.csv", stringsAsFactors = F)
cmap_experiments = read.csv("raw/cmap/cmap_drug_experiments.csv", stringsAsFactors =  F)

drug_preds = results[[1]] #from cmap_predictions.RData
dz_signature = results[[2]]

cmap_experiments_valid = merge(cmap_experiments, valid_instances, by="id")
#keep valid profiles; drugs (as long as matched to Drugbank)
cmap_experiments_valid = cmap_experiments_valid[cmap_experiments_valid$valid == 1 & cmap_experiments_valid$DrugBank.ID != "NULL", ]
 
drug_instances_all = merge(drug_preds, cmap_experiments_valid, by.x="exp_id", by.y="id")
write.csv(drug_instances_all, paste(cancer, "/drug/cmap_drug_predictions.csv", sep=""))

#sig cmap cutoff
cutoff = max(drug_instances_all$cmap_score[drug_instances_all$q < 0.01 & drug_instances_all$cmap_score < 0])

#visualize score distritbution
segmentWidth = 0.025
jitterWidth = 0.01
jitterHight = 0.01
color = ""
drug_instances_all$group = -1
drug_instances_all$score = (drug_instances_all$cmap_score)
pdf(paste(cancer, "/drug/cmap_final","_cmap_distribution.pdf", sep=""))
  createGGPlotDist(drug_instances_all, measureVar="score", groupVars = "group", main = "", ylab = "", xlab = "")
dev.off()

drug_instances = subset(drug_instances_all, q < 0.05)
drug_instances = drug_instances[order(drug_instances$cmap_score), ]

drug_instances_id <- drug_instances$exp_id[1:20] + 1 #the first column is the gene id
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

re_rank_instance = 1
if (re_rank_instance > 0){
  col_sorted = sort(cor(drug_dz_signature, method="spearman")["value",-1])    
  drug_dz_signature = drug_dz_signature[,c("value", names(col_sorted))]
}

drug_names <- sapply(2:ncol(drug_dz_signature), function(id){
  #need to minus 1 as in cmap_signatures, V1 is gene id.
  new_id = strtoi(paste(unlist(strsplit(as.character(colnames(drug_dz_signature)[id]),""))[-1], collapse="")) - 1 
  cmap_experiments_valid$name[cmap_experiments_valid$id == new_id]
})

in_vitro_list = c("pyrvinium", "niclosamide", "mebendazole", "bithionol")
pdf(paste(cancer, "/drug/", "cmap_predictions.pdf", sep=""))
  colPal <- redgreen(100)
  par(mar=c(12, 4, 2, 0.5))
  #image(t(druggable_targets_pos), col=redblue(2))
  image(t(drug_dz_signature), col= colPal,   axes=F, srt=45)
  axis(1,  at=seq(0,1,length.out= ncol( drug_dz_signature ) ), labels= F)
  cols = ifelse(tolower(drug_names) %in% in_vitro_list, "darkgreen", "black")
  text(x = seq(0,1,length.out=ncol( drug_dz_signature ) ), c(-0.05),
     labels = c( "HCC",drug_names), srt = 45, pos=2,offset=0.05, xpd = TRUE, cex=1.4, col = c("red", cols))
dev.off()
write.csv(drug_dz_signature, paste(cancer, "/drug/cmap_drug_dz_signature.csv", sep=""))

colnames(drug_dz_signature) = c( "HCC",drug_names)
#pheatmap((drug_dz_signature), col = colPal, cluster_row = FALSE, cluster_col = F, show_rownames = F, legend=T, filename = paste(cancer, "/drug/cmap_predictions_pheatmap.pdf", sep=""))
pheatmap(t(drug_dz_signature), col = colPal, cluster_row = FALSE, cluster_col = F, show_colnames = F, legend=F, filename = paste(cancer, "/drug/cmap_predictions_pheatmap_reverse.pdf", sep=""))

#cluster drugs 
drug_dz_signature_cluster = drug_dz_signature[, -c(1)]
colnames(drug_dz_signature_cluster) = drug_names
#pheatmap(drug_dz_signature_cluster, col = colPal, cluster_row = FALSE, clustering_distance_cols = "correlation", show_rownames = F, legend=F,  filename = paste(cancer, "/drug/cluster_by_drug.pdf", sep=""))

#cluster genes; require the signature includes gene symbol and gene id
gene_annot <- subset(dz_signature, select = c("Symbol", "GeneID"))
gene_ids_annot <- merge(gene_ids, gene_annot, by.x=1, by.y="GeneID", sort=F)
rownames(drug_dz_signature_cluster) = gene_ids_annot$Symbol
#pheatmap(drug_dz_signature_cluster, col = colPal, cellheight = 12, clustering_distance_rows = "correlation", show_rownames = T, legend=F,  filename = paste(cancer, "/drug/cluster_by_gene.pdf", sep=""))

#visualize reversed gene
drug_dz_signature_reversed = (drug_dz_signature - drug_dz_signature[, 1]) #abs
drug_dz_signature_reversed = drug_dz_signature_reversed[, -c(1)]
gene_names = as.character(gene_ids_annot$Symbol)

rownames(drug_dz_signature_reversed) = gene_names
colnames(drug_dz_signature_reversed) = drug_names
my.cols <- brewer.pal(9, "Blues")
#pheatmap(drug_dz_signature_reversed, col = my.cols , cellheight = 12,  show_rownames = T, legend=T,filename =  paste(cancer, "/drug/reversibility.pdf", sep="")) #,  

#transpose, let gene in the column
direction_matrix = drug_dz_signature - drug_dz_signature[, 1]
direction_matrix = sign(direction_matrix[, c(-1)])
annotation = data.frame(type=sign(as.vector(apply(direction_matrix, 1, sum))))
rownames(annotation) = gene_names
annotation$type = factor(annotation$type)
drug_dz_signature_rotate = t(drug_dz_signature_reversed)
pheatmap(-drug_dz_signature_rotate, col = redblue(100) , cellheight = 12,  cellwidth = 10, show_rownames = T, legend=T,  filename =  paste(cancer, "/drug/cmap_gene_reversed_two_color.pdf", sep="")) #,  filename = "reversibility.pdf"
pheatmap(abs(drug_dz_signature_rotate), col = my.cols , annotation = annotation, cellheight = 12,  cellwidth = 10, show_rownames = T, legend=T,  filename =  paste(cancer, "/drug/cmap_gene_reversed_one_color.pdf", sep="")) #,  filename = "reversibility.pdf"

##################
#violin plot 
drugs = drug_names #unique(drug_instances_all$name[drug_instances_all$q_value < 0.05 & drug_instances_all$cmap_score < 0 & tolower(drug_instances_all$name) %in% c(tolower(drugbank_drugs), "pyrvinium")])

#signicant cutoff
sig_cutoff = max(drug_instances_all$cmap_score[drug_instances_all$q < 0.05 & drug_instances_all$cmap_score < 0])

medians = NULL
means = NULL
ratios = NULL
counts = NULL
for (i in 1:length(drugs)){
  cmap_predictions_drug = subset(drug_instances_all, name == drugs[i])
  medians = c(medians, median(cmap_predictions_drug$cmap_score))
  means = c(means, mean(cmap_predictions_drug$cmap_score))
  ratios = c(ratios, nrow(subset(cmap_predictions_drug, cmap_score < sig_cutoff & q < 0.05)) / nrow(cmap_predictions_drug))
  counts = c(counts, nrow(cmap_predictions_drug))
}

drugs_rank = data.frame(drugs, medians, means, ratios, counts)
#rank by median
drugs_rank = drugs_rank[order(drugs_rank$medians, decreasing=F),]
#drugs_rank = subset(drugs_rank, counts > 1)
head(drugs_rank, 25)
tail(drugs_rank, 25)

#visualize top hits using violin plot
candidates = head(drugs_rank, 20)$drugs
hits_scores = subset(drug_instances_all, name %in% candidates, select = c("exp_id", "cmap_score", "name"))
names(hits_scores) = c("exp_id", "score", "group")
hits_scores$class = hits_scores$group
scores = rbind(hits_scores) #, all_scores
scores = merge(scores, subset(drug_instances_all, select=c("exp_id",  "cell_line")), by = "exp_id")
jitterWidth = 0.1
jitterHight = 0
color = ""
## data.frame must have at least three columns - score, group, class.
## In this case, data.frame$class = data.frame$group
pdf(paste(cancer, "/drug/cmap_prediction_violin.pdf", sep=""))
  drug_order <- aggregate(score ~ group + class, scores, median)
  drug_ordered <- drug_order$group[order(drug_order$score, decreasing=T)]
  scores$group = factor(scores$group, levels = c( drug_ordered))
  scores$cell_line = scores$cell_line
  createGGPlot(scores, measureVar="score", groupVars = "group", main = "", ylab = "Score in CMap", xlab = "", sig_cutoff)
dev.off()


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

drug_preds_sig = merge(drug_preds, lincs_experiments, by.x="exp_id", by.y="id")

lincs_predictions_all = subset(drug_preds_sig, cell_id %in% c("HEPG2", "HUH7") & DrugBank.ID != "NULL")
lincs_predictions_all = lincs_predictions_all[order(lincs_predictions_all$cmap_score),]

write.csv(lincs_predictions_all, paste(cancer, "/drug/lincs_drug_predictions.csv", sep=""))

#visualize score distritbution
segmentWidth = 0.025
jitterWidth = 0.01
jitterHight = 0.01
color = ""
lincs_predictions_all$group = -1
lincs_predictions_all$score = (lincs_predictions_all$cmap_score)
pdf(paste(cancer, "/drug/lincs_final","_lincs_distribution.pdf", sep=""))
  createGGPlotDist(lincs_predictions_all, measureVar="score", groupVars = "group", main = "", ylab = "", xlab = "")
dev.off()

lincs_predictions = lincs_predictions_all

#visualize top 20 drugs
#get drug signatures
sig_ids = lincs_predictions$exp_id[1:20]
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

in_vitro_list = c("pyrvinium", "niclosamide", "mebendazole", "bithionol")
###visualze top 20 genes
pdf(paste(cancer, "/drug/lincs_predictions.pdf", sep=""))
  layout(matrix(1))
  par(mar=c(12, 4, 1, 0.5))
  colPal <- redgreen(100)
  image(t(drug_dz_signature), col= colPal,   axes=F, srt=45)
  cols = ifelse(tolower(drug_names) %in% in_vitro_list, "darkgreen", "black")
  axis(1,  at=seq(0,1,length.out= ncol( drug_dz_signature ) ), labels= F)
  text(x = seq(0,1,length.out=ncol( drug_dz_signature ) ), c(-0.05),
       labels = c( "HCC",drug_names), srt = 45, pos=2,offset=0.05, xpd = TRUE, cex=1.4, col = c("red",   cols = ifelse(tolower(drug_names) %in% in_vitro_list, "darkgreen", "black")))
dev.off()
write.csv(drug_dz_signature, paste(cancer, "/drug/lincs_drug_dz_signature.csv", sep=""))

drug_dz_signature_reversed = (drug_dz_signature - drug_dz_signature[, 1]) #abs
drug_dz_signature_reversed = drug_dz_signature_reversed[, -c(1)]
gene_ids_annot <- merge(gene_ids, gene_annot, by.x=1, by.y="GeneID", sort=F)
gene_names = as.character(gene_ids_annot$Symbol)

rownames(drug_dz_signature_reversed) = as.character(gene_names)
colnames(drug_dz_signature_reversed) = drug_names
my.cols <- brewer.pal(9, "Blues")
pheatmap(drug_dz_signature_reversed, col = my.cols , cellheight = 12,  show_rownames = T, legend=T,filename = paste(cancer, "/drug/lincs_reversibility.pdf",sep="")) #,  

#transpose, let gene in the column
direction_matrix = drug_dz_signature - drug_dz_signature[, 1]
direction_matrix = sign(direction_matrix[, c(-1)])
annotation = data.frame(type=sign(as.vector(apply(direction_matrix, 1, sum))))
rownames(annotation) = gene_names
annotation$type = factor(annotation$type)

drug_dz_signature_rotate = t(drug_dz_signature_reversed)
pheatmap(-drug_dz_signature_rotate, col = redblue(100) , cellheight = 12,  cellwidth = 10, show_rownames = T, legend=T,  filename = paste(cancer, "/drug/lincs_gene_reversed_two_color.pdf", sep="")) #,  filename = "reversibility.pdf"
pheatmap(abs(drug_dz_signature_rotate), col = my.cols , annotation = annotation, cellheight = 12,  cellwidth = 10, show_rownames = T, legend=T,  filename = paste(cancer, "/drug/lincs_gene_reversed_one_color.pdf", sep="")) #,  filename = "reversibility.pdf"

###############
###enrichment analysis
drugs = unique(lincs_predictions$pert_iname[lincs_predictions$q < 0.05])
#signicant cutoff
sig_cutoff = max(lincs_predictions$cmap_score[lincs_predictions$q < 0.05 & lincs_predictions$cmap_score < 0])

medians = NULL
means = NULL
ratios = NULL
counts = NULL
for (i in 1:length(drugs)){
  lincs_predictions_drug = subset(lincs_predictions, pert_iname == drugs[i])
  medians = c(medians, median(lincs_predictions_drug$cmap_score))
  means = c(means, mean(lincs_predictions_drug$cmap_score))
  ratios = c(ratios, nrow(subset(lincs_predictions_drug, cmap_score < sig_cutoff & q < 0.05)) / nrow(lincs_predictions_drug))
  counts = c(counts, nrow(lincs_predictions_drug))
}

drugs_rank = data.frame(drugs, medians, means, ratios, counts)
drugs_rank = drugs_rank[order(drugs_rank$medians, decreasing=F),]
head(drugs_rank, 25)
tail(drugs_rank, 25)



#candidates are top 20 ranked drugs or the hits with top average score
drug_names = lincs_predictions$pert_iname[1:20]
candidates = drug_names #as.character(drugs_rank$drugs[1:20]) #c("sotalol", "mebendazole", "brimonidine", "epibatidine")
##################
#violin plot 
hits_scores = subset(lincs_predictions, pert_iname %in% candidates, select = c("exp_id", "cmap_score", "pert_iname"))
names(hits_scores) = c("exp_id", "score", "group")
hits_scores$class = hits_scores$group

scores = rbind(hits_scores) #, all_scores
scores = merge(scores, subset(lincs_predictions, select=c("exp_id",  "cell_id")), by = "exp_id")

segmentWidth = 0.025
jitterWidth = 0.1
jitterHight = 0
color = ""
## data.frame must have at least three columns - score, group, class.
## In this case, data.frame$class = data.frame$group
par(mar=c(4, 4, 2, 4))
pdf(paste(cancer, "/drug/lincs_prediction_violin.pdf", sep=""))
  drug_order <- aggregate(score ~ group + class, scores, median)
  drug_ordered <- drug_order$group[order(drug_order$score, decreasing=T)]
  scores$group = factor(scores$group, levels = drug_ordered)
  scores$cell_line = scores$cell_id
  createGGPlot(scores, measureVar="score", groupVars = "group", main = "", ylab = "Score in LINCS", xlab = "", sig_cutoff)
dev.off()


############
#merge LINCS and CMap
############

cmap_results = read.csv(paste(cancer, "/drug/cmap_drug_predictions.csv", sep=""), stringsAsFactors = F)
lincs_results = read.csv(paste(cancer, "/drug/lincs_drug_predictions.csv", sep=""), stringsAsFactors = F)

cmap_drugs = subset(cmap_results, q < 0.01 & cmap_score < 0)
lincs_drugs = subset(lincs_results, q < 0.01 & cmap_score < 0)

common_drugs = intersect(tolower(cmap_drugs$name), tolower(lincs_drugs$pert_iname))

write.csv(common_drugs, paste(cancer, "/drug/cmap_lincs_results_common.csv", sep=""))

source("http://faculty.ucr.edu/~tgirke/Documents/R_BioCond/My_R_Scripts/overLapper.R") # Imports required functions.
setlist <- list(CMAP = unique(tolower(cmap_drugs$name)), LINCS=unique(lincs_drugs$pert_iname)) # ,C=unique(set1$GeneID[set1$up_down=="down"]), D=unique(set2$GeneID[set2$up_down=="down"]))
setlist2 <- setlist[c(1,2)]; 
list2 <- overLapper(setlist=setlist2, sep="_", type="vennsets")
counts <- sapply(list2$Venn_List, length);
pdf(paste(cancer, "/drug/cmap_lincs_overlap.pdf", sep=""))
  vennPlot(counts=counts)
dev.off()

