subgroup = "alchohol"
setwd(paste("/Users/binchen1/Documents/stanford/hcc/data/LIHC/drug_subgroup/", subgroup, sep=""))

get.drug.name <- function(con, id, only_name=T, source="cmap"){
  if ( source == "cmap"){ #indicate it is from cmap
    id_new = strtoi(paste(unlist(strsplit(id,""))[-1], collapse="")) - 1
    query <- paste("select * from cmap_experiments where id =", id_new )
    rs <- dbSendQuery(con, query)        
    drug_info <- fetch(rs, n = 1)
    dbClearResult(rs)
    if (only_name){
      return(drug_info[1,"name"])
    }else{
      return(paste(drug_info[1,c("name","cell_line","duration","concentration")], collapse="_")) #)
    }
  }else{
    query = paste("select * from proj_lincs.sig_info where id ='", id, "'", sep="")
    rs <- dbSendQuery(con, query)        
    drug_info <- fetch(rs, n = 1)
    dbClearResult(rs)
    if (only_name){
      return(drug_info[1,"pert_iname"])
    }else{
      return(paste(drug_info[1,c("pert_iname","cell_id","pert_dose","pert_time")], collapse="_")) #
    }     
  }
}

get.landmark.info <- function(con){
  rs <- dbSendQuery(con, "select * from proj_lincs.probe_id_info where pool_id like '%epsilon%'")    
  probe_info <- fetch(rs, n = -1)
  dbClearResult(rs)
  return(probe_info)
}

get.geneid.from.synonym <- function(con, synonym){
  rs <- dbSendQuery(con, paste("select geneid from annot_gene.gene_synonym where tax_id = 9606 and synonym = '", synonym, "'", sep=""))    
  gene_info <- fetch(rs, n = -1)
  dbClearResult(rs)
  if (nrow(gene_info) > 0){
    return (gene_info[1,1])
  }else{
    return (NA)
  }
}

get.instance.sig <- function(id, con, version = "new", landmark=T){
  
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
  
  
  return (instance.sig)
}

createGGPlot <- function(currPheno, measureVar, groupVars, main="GGPLot", ylab="pass ylab value", xlab="pass xlab value") {
  temp = summarySE(currPheno, measurevar=measureVar, groupvars=groupVars)[, c(groupVars, measureVar, "se")]
  a = ggplot() +
    geom_violin(data=currPheno, aes(x=group,y=score), fill='grey',trim=F) +
    coord_flip() +
    geom_jitter(data=currPheno, aes(x=group,y=score, color=color), size = 2, position=position_jitter(width=jitterWidth, height=jitterHight)) +
    guides(color=FALSE) +
    xlab(xlab) +
    ylab(ylab) + 
    ggtitle(main) +
    theme(text = element_text(size=24)) +
    #theme(legend.position=c(0.5,0.9), legend.title=element_blank(), legend.text=element_text(size=24)) +
    theme(axis.text.x = element_text(size=24, angle = 90, hjust = 1), axis.text.y = element_text(size=24)) +
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
    geom_hline(yintercept = 0, linetype = "longdash")
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

##############################
######
library(pheatmap)
library(gplots)
library(ggplot2)
library("RMySQL")
mysql_drvr <-dbDriver("MySQL")
con <- dbConnect(mysql_drvr, group="client",host="buttelab-db1.stanford.edu",dbname="proj_repositioning")

signatures = read.table(paste("LIHC", subgroup, "_dz_signature_lincs.txt", sep=""), header=T,sep="\t") #at least cover gene id and value that can be either fold change or p value
signatures$value = - signatures$value #norm vs control

#############
#read drug-bank

drugbank_drugtargets = dbReadTable(con, "user_binchen1.v_drugbank_drug_targets_new")
#drugbank_drugtargets = subset(drugbank_drugtargets, gene_symbol != '', select = c("drug_name", "gene_symbol", "action"))
drugbank_drugtargets = subset(drugbank_drugtargets, type == "target")

query = "select dbid, GROUP_CONCAT(atc) as atcs from user_binchen1.drugbank_atcs_new group by dbid"
rs = dbSendQuery(con, query)
drug_atcs = fetch(rs, n=-1)

signatures_drugbank = merge(signatures, drugbank_drugtargets, by.x = "Symbol", by.y = "gene_symbol")
unique(signatures_drugbank$Gene)
write.table(signatures_drugbank, "LIHC_mapped_drugbank_new.txt", sep="\t", quote=F, col.names=T, row.names=F )


#############
###CMAP
library("gplots")
load('/Users/binchen1/Documents/stanford/lincs/data/cmap_signatures_updated.RData')

signatures = read.table(paste("LIHC", subgroup, "_dz_signature_cmap.txt", sep=""), header=T,sep="\t") #at least cover gene id and value that can be either fold change or p value
signatures$value = - signatures$value #norm vs control


query <- paste("select experiment_id, name from proj_repositioning.v_cmap_hits_valid  where valid=1 and subset_comparison_id = 'hcc_hbv' order by cmap_score  "    )         
rs <- dbSendQuery(con, query)    
drug_instances <- fetch(rs, n = -1)
dbClearResult(rs)

#only keep drugs
drugbank_drugs = dbReadTable(con, "user_binchen1.drugbank_drugs_new")
drugbank_drugs = drugbank_drugs$name
drug_instances = subset(drug_instances, tolower(name) %in% tolower(drugbank_drugs), select="experiment_id")

drug_instances_id <- drug_instances$experiment_id[1:20] + 1 #the first column is the gene id
#get candidate drugs
drug_signatures <- cmap_signatures[,c(1, drug_instances_id)] #the first column is the gene id

#load dz genes
dz_signature <- signatures #read.table("LIHC_dz_signature_cmap.txt",header=T,sep="\t") #at least cover gene id and value that can be either fold change or p value
dz_signature <- subset(dz_signature, !is.na(GeneID), select=c("GeneID", "value"))

drug_dz_signature <- merge(dz_signature, drug_signatures, by.x = "GeneID", by.y="V1")
drug_dz_signature <- drug_dz_signature[order(drug_dz_signature$value),]

#rerank 
drug_dz_signature[,2] = -drug_dz_signature[,2] # the higher rank indicate the more overexpressed, so we need to reverse 
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
  get.drug.name(con, colnames(drug_dz_signature)[id], T)
})


#drug_dz_signature = cbind(druggable_targets_pos, drug_dz_signature )

pdf("cmap_predictions.pdf")
  colPal <- redgreen(100)
  par(mar=c(12, 4, 2, 0.5))
  #image(t(druggable_targets_pos), col=redblue(2))
  image(t(drug_dz_signature), col= colPal,   axes=F, srt=45)
  axis(1,  at=seq(0,1,length.out= ncol( drug_dz_signature ) ), labels= F)
  text(x = seq(0,1,length.out=ncol( drug_dz_signature ) ), c(-0.05),
     labels = c( "HCC",drug_names), srt = 45, pos=2,offset=0.05, xpd = TRUE, cex=0.8)

dev.off()

pdf("cmap_predictions_cor_banner.pdf")
  layout(matrix(c(1,1,2,2), 2, 2, byrow = TRUE), heights= c(1,10))
  par(mar=c(0, 4, 2, 0.5))
  cors = cor(drug_dz_signature, method="spearman")["value",-1]
  image((as.matrix(c(NA, -cors))), col= brewer.pal(9, "Blues") ,   axes=F, srt=45) #brewer.pal(length(cors), "Blues")
  
  par(mar=c(12, 4, 0, 0.5))
  colPal <- redgreen(100)
  image(t(drug_dz_signature), col= colPal,   axes=F, srt=45)
  axis(1,  at=seq(0,1,length.out= ncol( drug_dz_signature ) ), labels= F)
  text(x = seq(0,1,length.out=ncol( drug_dz_signature ) ), c(-0.05),
     labels = c( "HCC",drug_names), srt = 45, pos=2,offset=0.05, xpd = TRUE, cex=0.8)
dev.off()

colnames(drug_dz_signature) = c( "HCC",drug_names)
pheatmap((drug_dz_signature), col = colPal, cluster_row = FALSE, cluster_col = F, show_rownames = F, legend=T, filename = "cmap_predictions_pheatmap.pdf")
pheatmap(t(drug_dz_signature), col = colPal, cluster_row = FALSE, cluster_col = F, show_colnames = F, legend=F, filename = "cmap_predictions_pheatmap_reverse.pdf")


#cluster drugs 
drug_dz_signature_cluster = drug_dz_signature[, -c(1)]
colnames(drug_dz_signature_cluster) = drug_names
pheatmap(drug_dz_signature_cluster, col = colPal, cluster_row = FALSE, clustering_distance_cols = "correlation", show_rownames = F, legend=F)

#cluster genes
gene_annot <- signatures #read.table("LIHC_dz_signature_cmap.txt",header=T,sep="\t") 
gene_annot <- subset(gene_annot, select = c("Symbol", "GeneID"))

gene_ids_annot <- merge(gene_ids, gene_annot, by.x=1, by.y="GeneID", sort=F)
rownames(drug_dz_signature_cluster) = gene_ids_annot$Symbol
pheatmap(drug_dz_signature_cluster, col = colPal, cellheight = 12, clustering_distance_rows = "correlation", show_rownames = T, legend=F,  filename = "cluster_by_gene.pdf")

#visualize reversed gene
mapped_targets = read.csv("LIHC_mapped_drugbank_new.txt",sep="\t")
targets = as.character(unique(mapped_targets$Symbol))  # c("GART", "TNFSF13B", "PDE4D", "THRA", "VAMP2", "SMOX")
drug_dz_signature_reversed = (drug_dz_signature - drug_dz_signature[, 1]) #abs
drug_dz_signature_reversed = drug_dz_signature_reversed[, -c(1)]
gene_names = as.character(gene_ids_annot$Symbol)
gene_names = sapply(gene_names, function(name){
  if (name %in% targets){
    paste(name, "*", sep="")
  }else{
    name
  }
})

rownames(drug_dz_signature_reversed) = gene_names
colnames(drug_dz_signature_reversed) = drug_names
library(RColorBrewer)
my.cols <- brewer.pal(9, "Blues")
pheatmap(drug_dz_signature_reversed, col = my.cols , cellheight = 12,  show_rownames = T, legend=T,filename = "reversibility.pdf") #,  

#transpose, let gene in the column
direction_matrix = drug_dz_signature - drug_dz_signature[, 1]
direction_matrix = sign(direction_matrix[, c(-1)])
annotation = data.frame(type=sign(as.vector(apply(direction_matrix, 1, sum))))
rownames(annotation) = gene_names
annotation$type = factor(annotation$type)

drug_dz_signature_rotate = t(drug_dz_signature_reversed)
pheatmap(-drug_dz_signature_rotate, col = redblue(100) , cellheight = 12,  cellwidth = 10, show_rownames = T, legend=T,  filename = "gene_reversed_two_color.pdf") #,  filename = "reversibility.pdf"
pheatmap(abs(drug_dz_signature_rotate), col = my.cols , annotation = annotation, cellheight = 12,  cellwidth = 10, show_rownames = T, legend=T,  filename = "gene_reversed_one_color.pdf") #,  filename = "reversibility.pdf"

#redblue(100)

#reverse targets
drug_dz_signature_reversed = drug_dz_signature - drug_dz_signature[, 1]
drug_dz_signature_reversed = cbind(gene_ids, drug_dz_signature_reversed)
drug_dz_signature_reversed_target = drug_dz_signature_reversed[ gene_ids %in% mapped_targets$GeneID, ]

drug_dz_signature_reversed_target = as.matrix(drug_dz_signature_reversed_target)
rownames(drug_dz_signature_reversed_target) = drug_dz_signature_reversed_target[,1]
drug_dz_signature_reversed_target = abs(drug_dz_signature_reversed_target[,-c(1,2)])
#image(t(drug_dz_signature_reversed_target), col= colPal,   axes=F, srt=45)

#plot reversed targets
pheatmap(drug_dz_signature_reversed_target, col = colPal, cellheight = 12, clustering_distance_rows = "correlation", show_rownames = T, legend=F,  filename = "cluster_by_gene.pdf")



####
##lincs
#lincs results 
#lincs results 
library("RMySQL")
mysql_drvr <-dbDriver("MySQL")
con <- dbConnect(mysql_drvr, group="client",host="buttelab-db1.stanford.edu",dbname="proj_lincs")

#genes shared by landmark genes and dz signatures
landmark = get.landmark.info(con)
signatures = read.table(paste("LIHC", subgroup, "_dz_signature_lincs.txt", sep=""), header=T,sep="\t") #at least cover gene id and value that can be either fold change or p value
signatures$value = - signatures$value #norm vs control

common_genes = merge(landmark, signatures, by.y="GeneID", by.x="gene_id")
write.csv(common_genes, "sig_metah_overlap_lincs.csv", row.names=F)
common_genes = common_genes[order(common_genes$value), ]
dz_sig = subset(common_genes, select=c("gene_id", "value"))
#cell lines related to dz

rs <- dbSendQuery(con, paste("select * from proj_lincs.v_lincs_predictions_drugs where cell_id in ('HEPG2', 'HUH7') and subset_comparison_id ='", "hcc_hbv", "'", sep=""))    
lincs_predictions <- fetch(rs, n = -1)
lincs_predictions = lincs_predictions[order(lincs_predictions$cmap_score),]

#visualize top 20 drugs
#get drug signatures
sig_ids = lincs_predictions$id[1:20]
drug_names = lincs_predictions$pert_iname[1:20]
sigs <- data.frame(gene_id =landmark$gene_id)
for (sig_id in sig_ids){
  sig = get.instance.sig(sig_id, con)[,3]
  sigs[,as.character(sig_id)] = sig
}

drug_dz_signature = merge(dz_sig, sigs, by="gene_id")
drug_dz_signature = drug_dz_signature[order(drug_dz_signature$value, decreasing = T),]
gene_ids = drug_dz_signature[,1]
drug_dz_signature = drug_dz_signature[, -c(1)]

#since the value is z score. the higher z score should have a higher rank...so we reverse the score
drug_dz_signature = -1 * drug_dz_signature
for (i in 1:ncol(drug_dz_signature)){
  drug_dz_signature[,i] = rank(drug_dz_signature[,i])
}

re_rank_instance = 1
if (re_rank_instance > 0){
  col_sorted = sort(cor(drug_dz_signature, method="spearman")["value",-1])    
  drug_dz_signature = drug_dz_signature[,c("value", names(col_sorted))]
}

drug_names <- sapply(2:ncol(drug_dz_signature), function(id){
  get.drug.name(con, colnames(drug_dz_signature)[id], T, "lincs")
})


gene_annot <- signatures #read.table("LIHC_dz_signature_lincs.txt",header=T,sep="\t") 
gene_annot <- subset(gene_annot, select = c("Symbol", "GeneID"))


###visualze top 20 genes
pdf("lincs_predictions.pdf")
  layout(matrix(1))
  par(mar=c(12, 4, 1, 0.5))
  colPal <- redgreen(100)
  image(t(drug_dz_signature), col= colPal,   axes=F, srt=45)
  axis(1,  at=seq(0,1,length.out= ncol( drug_dz_signature ) ), labels= F)
  text(x = seq(0,1,length.out=ncol( drug_dz_signature ) ), c(-0.05),
       labels = c( "HCC",drug_names), srt = 45, pos=2,offset=0.05, xpd = TRUE, cex=0.8)

dev.off()

pdf("lincs_predictions_cor_banner.pdf")
  layout(matrix(c(1,1,2,2), 2, 2, byrow = TRUE), heights= c(1,10))
  par(mar=c(0, 4, 2, 0.5))
  cors = cor(drug_dz_signature, method="spearman")["value",-1]
  image((as.matrix(c(NA, -cors))), col= brewer.pal(length(cors), "Blues") ,   axes=F, srt=45) #brewer.pal(length(cors), "Blues")
  
  par(mar=c(12, 4, 0, 0.5))
  colPal <- redgreen(100)
  image(t(drug_dz_signature), col= colPal,   axes=F, srt=45)
  axis(1,  at=seq(0,1,length.out= ncol( drug_dz_signature ) ), labels=c( "HCC",drug_names),srt=45, las =2)
dev.off()

mapped_targets = read.csv("LIHC_mapped_drugbank_new.txt",sep="\t")
targets = as.character(unique(mapped_targets$Symbol))  # c("GART", "TNFSF13B", "PDE4D", "THRA", "VAMP2", "SMOX")
drug_dz_signature_reversed = (drug_dz_signature - drug_dz_signature[, 1]) #abs
drug_dz_signature_reversed = drug_dz_signature_reversed[, -c(1)]
gene_ids_annot <- merge(gene_ids, gene_annot, by.x=1, by.y="GeneID", sort=F)
gene_names = as.character(gene_ids_annot$Symbol)
gene_names = sapply(gene_names, function(name){
  if (name %in% targets){
    paste(name, "*", sep="")
  }else{
    name
  }
})

rownames(drug_dz_signature_reversed) = as.character(gene_names)
colnames(drug_dz_signature_reversed) = drug_names
library(RColorBrewer)
my.cols <- brewer.pal(9, "Blues")
pheatmap(drug_dz_signature_reversed, col = my.cols , cellheight = 12,  show_rownames = T, legend=T,filename = "lincs_reversibility.pdf") #,  

#transpose, let gene in the column
direction_matrix = drug_dz_signature - drug_dz_signature[, 1]
direction_matrix = sign(direction_matrix[, c(-1)])
annotation = data.frame(type=sign(as.vector(apply(direction_matrix, 1, sum))))
rownames(annotation) = gene_names
annotation$type = factor(annotation$type)

drug_dz_signature_rotate = t(drug_dz_signature_reversed)
pheatmap(-drug_dz_signature_rotate, col = redblue(100) , cellheight = 12,  cellwidth = 10, show_rownames = T, legend=T,  filename = "lincs_gene_reversed_two_color.pdf") #,  filename = "reversibility.pdf"
pheatmap(abs(drug_dz_signature_rotate), col = my.cols , annotation = annotation, cellheight = 12,  cellwidth = 10, show_rownames = T, legend=T,  filename = "lincs_gene_reversed_one_color.pdf") #,  filename = "reversibility.pdf"

###############
###enrichment analysis

drug_names = lincs_predictions$pert_iname[1:20]

drugs = unique(lincs_predictions$pert_iname[lincs_predictions$p_value < 0.05])

#signicant cutoff
sig_cutoff = max(lincs_predictions$cmap_score[lincs_predictions$p_value < 0.05 & lincs_predictions$cmap_score < 0])

medians = NULL
means = NULL
ratios = NULL
counts = NULL
for (i in 1:length(drugs)){
  lincs_predictions_drug = subset(lincs_predictions, pert_iname == drugs[i])
  #lines(density(lincs_predictions_drug$cmap_score), col=i)
  medians = c(medians, median(lincs_predictions_drug$cmap_score))
  means = c(means, mean(lincs_predictions_drug$cmap_score))
  ratios = c(ratios, nrow(subset(lincs_predictions_drug, cmap_score < sig_cutoff & p_value < 0.05)) / nrow(lincs_predictions_drug))
  counts = c(counts, nrow(lincs_predictions_drug))
}

drugs_rank = data.frame(drugs, medians, means, ratios, counts)
drugs_rank = drugs_rank[order(drugs_rank$medians, decreasing=F),]
drugs_rank = subset(drugs_rank, counts > 1)
head(drugs_rank, 25)
tail(drugs_rank, 25)

library(plyr)
cdata <- ddply(lincs_predictions, .(pert_iname), summarise, 
               mean_score = mean(cmap_score),
               len = length(cmap_score)
)
qnorm(0.05, mean = mean(cdata$mean_score), sd = sd(cdata$mean_score), lower.tail = T )


par(mar=c(4, 4, 2, 4))

pdf("lincs_predictions_density.pdf")
plot(density(lincs_predictions$cmap_score), col = 1, xlab="cmap score", main="")
candidates = as.character(drugs_rank$drugs[1:10]) #c("sotalol", "mebendazole", "brimonidine", "epibatidine")
for (i in 1:length(candidates)){
  lincs_predictions_drug = subset(lincs_predictions, pert_iname == candidates[i])
  lines(density(lincs_predictions_drug$cmap_score), col= i + 1)
}
legend("topright", c("all drugs", candidates), col = c(1:(length(candidates)+1)), lty=1)
dev.off()
#pheatmap(as.matrix(lincs_predictions$cmap_score), cluster_row = F, cluster_col = F, col = redgreen(1000),show_rownames = T, legend=T)




##################
#violin plot 
hits_scores = subset(lincs_predictions, pert_iname %in% candidates, select = c("cmap_score", "pert_iname"))
names(hits_scores) = c("score", "group")
hits_scores$class = hits_scores$group

all_scores = subset(lincs_predictions,  select = c("cmap_score"))
names(all_scores) = c("score")
all_scores$group = "ALL"
all_scores$class = all_scores$group

scores = rbind(hits_scores, all_scores)

segmentWidth = 0.025
jitterWidth = 0.1
jitterHight = 0
color = ""
library("ggplot2")
library("plyr")
## data.frame must have at least three columns - score, group, class.
## In this case, data.frame$class = data.frame$group


pdf("lincs_prediction_violin.pdf")
drug_order <- aggregate(score ~ group + class, scores, mean)
drug_ordered <- drug_order$group[order(drug_order$score, decreasing=T)]
scores$group = factor(scores$group, levels = drug_ordered)

createGGPlot(scores, measureVar="score", groupVars = "group", main = "", ylab = "cmap score", xlab = "")
dev.off()


#find common drugs 
query <- paste("select * from proj_repositioning.v_cmap_hits_valid  where cmap_score < 0 and valid=1 and subset_comparison_id = 'LIHC_cmap' and q_value < 0.001 order by cmap_score  "    )         
rs <- dbSendQuery(con, query)    
cmap_results <- fetch(rs, n = -1)
drugbank_drugs = dbReadTable(con, "user_binchen1.drugbank_drugs_new")
drugbank_drugs = drugbank_drugs$name
cmap_results = subset(cmap_results, tolower(name) %in% tolower(drugbank_drugs))
write.csv(cmap_results, "cmap_results.csv")

rs <- dbSendQuery(con, paste("select experiment_id, cmap_score, q_value,q_value, sig_id,pert_desc, pert_time,pert_dose, pert_iname, cell_id from proj_lincs.v_lincs_predictions_drugs where cmap_score<0 and q_value < 0.001 and subset_comparison_id ='", "LIHC_lincs_landmark", "' order by cmap_score", sep=""))    
lincs_results <- fetch(rs, n = -1)
write.csv(lincs_results, "lincs_results.csv")

rs <- dbSendQuery(con, paste("select experiment_id, cmap_score, q_value,q_value, sig_id,pert_desc, pert_time,pert_dose, pert_iname, cell_id from proj_lincs.v_lincs_predictions_all where cell_id in ('HEPG2', 'HUH7') and cmap_score<0 and q_value < 0.001 and subset_comparison_id ='", "LIHC_lincs_landmark", "' order by cmap_score", sep=""))    
lincs_results_HCC_cell_line <- fetch(rs, n = -1)
write.csv(lincs_results_HCC_cell_line, "lincs_results_HCC_cell_line.csv")

common_drugs = intersect(tolower(cmap_results$name), tolower(lincs_results_HCC_cell_line$pert_iname))

write.csv(common_drugs, "cmap_lincs_results_HCC_cell_line_common.csv")

source("http://faculty.ucr.edu/~tgirke/Documents/R_BioCond/My_R_Scripts/overLapper.R") # Imports required functions.
setlist <- list(CMAP = unique(tolower(cmap_results$name)), LINCS=unique(lincs_results_HCC_cell_line$pert_iname)) # ,C=unique(set1$GeneID[set1$up_down=="down"]), D=unique(set2$GeneID[set2$up_down=="down"]))
setlist2 <- setlist[c(1,2)]; 
list2 <- overLapper(setlist=setlist2, sep="_", type="vennsets")
counts <- sapply(list2$Venn_List, length);
pdf("cmap_lincs_overlap.pdf")
vennPlot(counts=counts)
dev.off()

