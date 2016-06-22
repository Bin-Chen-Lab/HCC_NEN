setwd("/Users/binchen1/Documents/stanford/hcc/data/LIHC/drug/")

###function

cmap_score <- function(sig_up, sig_down, drug_signature) {
  #the old function does not support the input list with either all up genes or all down genes, this new function attempts to addess this, not fully validated
  #Note. I think that creating the anonymous functions in each iteration of the sapply's below is slowing things down. Predefine them eventually.
  num_genes = nrow(drug_signature)
  ks_up <- 0
  ks_down <- 0
  connectivity_score <- 0
  
  # I think we are re-ranking because the GeneID mapping changed the original rank range
  drug_signature[,"rank"] <- rank(drug_signature[,"rank"])
  
  # Merge the drug signature with the disease signature by GeneID. This becomes the V(j) from the algorithm description
  up_tags_rank <- merge(drug_signature, sig_up, by.x = "ids", by.y = 1)
  down_tags_rank <- merge(drug_signature, sig_down, by.x = "ids", by.y = 1)
  
  up_tags_position <- sort(up_tags_rank$rank)
  down_tags_position <- sort(down_tags_rank$rank)
  
  num_tags_up <- length(up_tags_position)
  num_tags_down <- length(down_tags_position)
  
  # 
  if(num_tags_up > 1) {
    a_up <- 0
    b_up <- 0
    
    a_up <- max(sapply(1:num_tags_up,function(j) {
      j/num_tags_up - up_tags_position[j]/num_genes
    }))
    b_up <- max(sapply(1:num_tags_up,function(j) {
      up_tags_position[j]/num_genes - (j-1)/num_tags_up
    }))
    
    if(a_up > b_up) {
      ks_up <- a_up
    } else {
      ks_up <- -b_up
    }
  }else{
    ks_up <- 0
  }
  
  if (num_tags_down > 1){
    
    a_down <- 0
    b_down <- 0
    
    a_down <- max(sapply(1:num_tags_down,function(j) {
      j/num_tags_down - down_tags_position[j]/num_genes
    }))
    b_down <- max(sapply(1:num_tags_down,function(j) {
      down_tags_position[j]/num_genes - (j-1)/num_tags_down
    }))
    
    if(a_down > b_down) {
      ks_down <- a_down
    } else {
      ks_down <- -b_down
    }
  }else{
    ks_down <- 0
  }
  
  if (ks_up == 0 & ks_down != 0){ #only down gene inputed
    connectivity_score <- -ks_down
  }else if (ks_up !=0 & ks_down == 0){ #only up gene inputed
    connectivity_score <- ks_up
  }else if (sum(sign(c(ks_down,ks_up))) == 0) {
    connectivity_score <- ks_up - ks_down # different signs
  }
  
  return(connectivity_score)
}

###########

library(gplots)
library(ggplot2)
require(XLConnect)
library(pheatmap)

source("http://faculty.ucr.edu/~tgirke/Documents/R_BioCond/My_R_Scripts/overLapper.R") # Imports required functions.


signatures <- read.table("LIHC_dz_signature_lincs_final.txt",header=T,sep="\t") #at least cover gene id and value that can be either fold change or p value
signatures2 <- read.table("LIHC_dz_signature_cmap_final_v1.txt",header=T,sep="\t") #at least cover gene id and value that can be either fold change or p value

#cluster drug expression usign disease signatures
data2 = read.csv("../../microarray//Group3/data2.csv")
data2_expr = subset(data2, GENE_ID %in% signatures$GeneID & fail.count<3, select=c("GeneSymbol", "N_576", "N_577", "N_578", "N_596", "N_597", "N_598"))
data2_expr = aggregate(. ~ GeneSymbol, data2_expr, mean)
rownames(data2_expr) = data2_expr$GeneSymbol
data2_expr = data2_expr[,-1]
annotation = data.frame(type=c("control", "control", "control", "NEN", "NEN", "NEN"))
rownames(annotation) = c("N_576", "N_577", "N_578" , "N_596", "N_597" , "N_598")
annotation$type = as.factor(annotation$type)
annotation = subset(annotation, select=c("type"))

my.cols <- greenred(100) #brewer.pal(9, "Blues")
pheatmap(scale(data2_expr), col = my.cols, annotation = annotation,   
         show_colnames=T, legend=T, show_rownames=F
)
#596 does not colaborate with 597 and 598

sig_data = merge(data.frame(Symbol = rownames(data2_expr), data2_expr), signatures, by = "Symbol")
cor(sig_data[,c("N_576", "N_577", "N_578" , "N_596", "N_597" , "N_598", "value")], method="spearman")

#########
#reverse disease gene expression
#pdx
data2 = read.csv("../../microarray//Group3/data2.csv")
data2 = subset(data2, fail.count>=0) #consider all counts, was <3
data2_gene = aggregate(Test.Control.fc ~ GENE_ID, data2, mean)
sig_data = merge(data2_gene, signatures, by.x="GENE_ID", by.y="GeneID")
cor.test(sig_data$Test.Control.fc, sig_data$value, method="spearman")

plot(sig_data$Test.Control.fc, sig_data$value, xlab = "NEN-mediated expression change", ylab="HCC expression")
#write.csv(sig_data, "../doc/niclosamide_manuscript/data/NEN_HCC_expression.csv")

#sig genes
data2_sig = subset(data2, fail.count<3 & abs(Test.Control.fc) >1.5 & Test.Control.adj.pval < 0.05)

data2_sig_dz = merge(data2_sig, signatures, by.x="GENE_ID", by.y="GeneID")

#hepg2
g2_data2 = read.csv("../../microarray//Group1/data2.csv")
g2_data2 = subset(g2_data2, fail.count>=0) #<2
g2_data2_gene_NEN = aggregate(G2.NEN.G2.fc  ~ GENE_ID, g2_data2, mean)
g2_data2_gene_NICLO = aggregate(G2.Niclo.G2.fc  ~ GENE_ID, g2_data2, mean)

sig_g2_NEN = merge(g2_data2_gene_NEN, signatures,  by.x="GENE_ID", by.y="GeneID")
cor.test(sig_g2_NEN$G2.NEN.G2.fc, sig_g2_NEN$value, method="spearman")

sig_g2_NICLO = merge(g2_data2_gene_NICLO, signatures,  by.x="GENE_ID", by.y="GeneID")
cor.test(sig_g2_NICLO$G2.Niclo.G2.fc, sig_g2_NICLO$value, method="spearman")


#huh7
H7_data2 = read.csv("../../microarray//Group2/data2.csv")
H7_data2 = subset(H7_data2, fail.count>=0) #>=0
H7_data2_gene_NEN = aggregate(H7.NEN.H7.fc  ~ GENE_ID, H7_data2, mean)
H7_data2_gene_NICLO = aggregate(H7.Niclo.H7.fc  ~ GENE_ID, H7_data2, mean)

sig_H7_NEN = merge(H7_data2_gene_NEN, signatures,  by.x="GENE_ID", by.y="GeneID")
cor.test(sig_H7_NEN$H7.NEN.H7.fc, sig_H7_NEN$value, method="spearman")

sig_H7_NICLO = merge(H7_data2_gene_NICLO, signatures,  by.x="GENE_ID", by.y="GeneID")
cor.test(sig_g2_NICLO$G2.Niclo.G2.fc, sig_g2_NICLO$value, method="spearman")

#hep40
Hep40_data2 = read.csv("../../microarray//Group4/data2.csv")
#Hep40_data2 = subset(Hep40_data2, fail.count>=0)
Hep40_data2_gene_NEN = aggregate(Hep40.NEN.Hep40.fc   ~ ENTREZ_GENE_ID, Hep40_data2, mean)
Hep40_data2_gene_NICLO = aggregate(Hep40.Niclo.Hep40.fc   ~ ENTREZ_GENE_ID, Hep40_data2, mean)
Hep40_data2_gene_NEN$GENE_ID = Hep40_data2_gene_NEN$ENTREZ_GENE_ID
Hep40_data2_gene_NICLO$GENE_ID = Hep40_data2_gene_NICLO$ENTREZ_GENE_ID

sig_Hep40_NEN = merge(Hep40_data2_gene_NEN, signatures,  by.x="GENE_ID", by.y="GeneID")
cor.test(sig_Hep40_NEN$Hep40.NEN.Hep40.fc, sig_Hep40_NEN$value, method="spearman")

sig_Hep40_NICLO = merge(Hep40_data2_gene_NICLO, signatures,  by.x="GENE_ID", by.y="GeneID")
cor.test(sig_Hep40_NICLO$Hep40.Niclo.Hep40.fc, sig_Hep40_NICLO$value, method="spearman")



in_vitro_all = merge(signatures,  g2_data2_gene_NEN, by.y="GENE_ID", by.x="GeneID", all.x=T)
in_vitro_all = merge(in_vitro_all, g2_data2_gene_NICLO,  by.y="GENE_ID", by.x="GeneID", all.x=T)
in_vitro_all = merge(in_vitro_all, H7_data2_gene_NEN,  by.y="GENE_ID", by.x="GeneID", all.x=T)
in_vitro_all = merge(in_vitro_all, H7_data2_gene_NICLO,  by.y="GENE_ID", by.x="GeneID", all.x=T)
in_vitro_all = merge(in_vitro_all, data2_gene,  by.y="GENE_ID", by.x="GeneID")
#in_vitro_all$PDX.invo = in_vitro_all$Test.Control.fc

cor.test(in_vitro_all$G2.NEN.G2.fc, in_vitro_all$G2.Niclo.G2.fc, method="spearman")
cor.test(in_vitro_all$H7.NEN.H7.fc, in_vitro_all$H7.Niclo.H7.fc, method = "spearman")


drug_dz_signature = subset(in_vitro_all, select=c("value",  "G2.Niclo.G2.fc", "G2.NEN.G2.fc", "H7.NEN.H7.fc", "H7.Niclo.H7.fc")) #"H7.NEN.H7.fc", "H7.Niclo.H7.fc", "Test.Control.fc"
drug_dz_signature_rank = drug_dz_signature
for (i in 1:ncol(drug_dz_signature)){
  drug_dz_signature_rank[, i] = rank(drug_dz_signature[,i])
}

drug_dz_signature_rank = drug_dz_signature_rank[order(drug_dz_signature_rank$value), ]

pdf("hcc_nen_niclosamide_reduced.pdf")
  colPal <- redgreen(100)
  par(mar=c(12, 4, 2, 0.5))
  #image(t(druggable_targets_pos), col=redblue(2))
  image(t(drug_dz_signature_rank), col= colPal,   axes=F, srt=45)
  axis(1,  at=seq(0,1,length.out= ncol( drug_dz_signature_rank ) ), labels= F)
  text(x = seq(0,1,length.out=ncol( drug_dz_signature_rank ) ), c(-0.05),
       labels = c( "HCC", "Niclosamide", "NEN"), srt = 45, pos=2,offset=0.05, xpd = TRUE, cex=1.4)
dev.off()

drug_dz_signature = subset(in_vitro_all, select=c("value",  "Test.Control.fc")) #"H7.NEN.H7.fc", "H7.Niclo.H7.fc", 
drug_dz_signature_rank = drug_dz_signature
for (i in 1:ncol(drug_dz_signature)){
  drug_dz_signature_rank[, i] = rank(drug_dz_signature[,i])
}

drug_dz_signature_rank = drug_dz_signature_rank[order(drug_dz_signature_rank$value), ]
pdf("hcc_nen_niclosamide_pdx_reduced.pdf")
  colPal <- redgreen(100)
  par(mar=c(12, 4, 2, 0.5))
  #image(t(druggable_targets_pos), col=redblue(2))
  image(t(drug_dz_signature_rank), col= colPal,   axes=F, srt=45)
  axis(1,  at=seq(0,1,length.out= ncol( drug_dz_signature_rank ) ), labels= F)
  text(x = seq(0,1,length.out=ncol( drug_dz_signature_rank ) ), c(-0.05),
       labels = c( "HCC", "NEN in vivo"), srt = 45, pos=2,offset=0.05, xpd = TRUE, cex=1.4)
dev.off()

#compute significance using GSEA
drug_gene_rank = data2_gene #H7_data2_gene_NEN H7_data2_gene_NICLO data2_gene
drug_gene_rank$rank = rank(-1 * drug_gene_rank$Test.Control.fc)
drug_gene_rank$ids = drug_gene_rank$GENE_ID

up = signatures$GeneID[signatures$value > 0]
down = signatures$GeneID[signatures$value < 0]
score = cmap_score(up, down, drug_gene_rank)

randoms = sapply(1:10000, function(id){
  #random select up/down genes
  drug_gene_rank$rank = sample(drug_gene_rank$rank, length(drug_gene_rank$rank))
  cmap_score(up, down, drug_gene_rank)  
})

q = sum(randoms < score)/length(randoms) #H7_data2_gene_NEN: 0.25,, 
q

######################
#find reversed genes
#fist find the rank of each gene in samples, and then random shuffle genes to get fake drug signatures
drug_dz_signature = subset(in_vitro_all, select=c("value",  "G2.Niclo.G2.fc", "G2.NEN.G2.fc","Test.Control.fc")) #, "Test.Control.fc": PDX data;  "H7.NEN.H7.fc", "H7.Niclo.H7.fc", 
drug_dz_signature_rank = drug_dz_signature
for (i in 1:ncol(drug_dz_signature)){
  drug_dz_signature_rank[, i] = rank(drug_dz_signature[,i])
}

real = drug_dz_signature_rank - drug_dz_signature_rank[,1]
real = real[, -1] #remove the first column
real_mean = apply(real, 1, mean)

permutations = 100000
permutation_m = matrix(rep(0, length(real_mean) * permutations), nrow=length(real_mean), ncol = permutations)
for (i in 1:permutations){
  print(i)
  #shuffle disease rank
  drug_dz_signature_rank_random = drug_dz_signature_rank
#  drug_dz_signature_rank_random[, 1] = sample(drug_dz_signature_rank_random[, 1], nrow(drug_dz_signature_rank_random))
  #shuffle drugs
  for (j in 2:ncol(drug_dz_signature_rank_random)){
    drug_dz_signature_rank_random[,j] = sample(drug_dz_signature_rank_random[,j], nrow(drug_dz_signature_rank_random))
  }
  random = drug_dz_signature_rank_random - drug_dz_signature_rank_random[,1]
  random = random[, -1] #remove the first column
  random_mean = apply(random, 1, mean)
  permutation_m[,i] = random_mean
}

p_up = sapply(1:length(real_mean), function(id){
  if (sum(permutation_m[id, ] > real_mean[id]) == 0){
    1/permutations #cannot be 0
  }else{
    sum(permutation_m[id, ] > real_mean[id])/permutations
  }
})
p_up_adj = p.adjust(p_up, "fdr")

p_down = sapply(1:length(real_mean), function(id){
  if (sum(permutation_m[id, ] < real_mean[id]) == 0){
    1/permutations #cannot be 0
  }else{
   sum(permutation_m[id, ] < real_mean[id])/permutations
  }
})
p_down_adj = p.adjust(p_down, "fdr")

gene_reversed = data.frame(GeneID = in_vitro_all$GeneID, Symbol = in_vitro_all$Symbol, p_up, p_up_adj, p_down, p_down_adj)
subset(gene_reversed, p_up < 0.05)
subset(gene_reversed, p_down < 0.05)
write.csv(gene_reversed, "../../microarray/reversed_genes_reduced.csv")


library(ggplot2)

gene_reversed_up = subset(gene_reversed, p_up_adj<0.25)
gene_reversed_up$log_q = log(gene_reversed_up$p_up_adj)
gene_reversed_up = gene_reversed_up[order(gene_reversed_up$p_up_adj),]
pdf("dz_induced_by_NEN_niclo_reduced.pdf")
ggplot(gene_reversed_up, aes(x=(Symbol), y = log_q)) +
  geom_bar(stat = "identity", fill="red") + 
  coord_flip() +
  theme(panel.background = element_rect(fill = 'white', colour = 'white'), axis.text = element_text(colour = "black", size=10)) +
  xlab("") +
  ylab("log(p value)") + 
  scale_fill_manual(values=c( "red", "green")) +
  theme(axis.text.x = element_text(size=17), axis.text.y = element_text(size=17), axis.title.y = element_text(size=18),  axis.title.x = element_text(size=18)) 
dev.off()

gene_reversed_down = subset(gene_reversed, p_down_adj<0.25)
gene_reversed_down$log_q = log(gene_reversed_down$p_down_adj)
gene_reversed_down = gene_reversed_down[order(gene_reversed_down$p_down_adj),]
pdf("dz_suppressed_by_NEN_niclo_reduced.pdf")
ggplot(gene_reversed_down, aes(x=(Symbol), y = log_q)) +
  geom_bar(stat = "identity", fill="green") + 
  coord_flip() +
  theme(panel.background = element_rect(fill = 'white', colour = 'white'), axis.text = element_text(colour = "black", size=10)) +
  xlab("") +
  ylab("log(p value)") + 
  scale_fill_manual(values=c( "red", "green")) +
  theme(axis.text.x = element_text(size=17), axis.text.y = element_text(size=17), axis.title.y = element_text(size=18),  axis.title.x = element_text(size=18)) 
dev.off()




