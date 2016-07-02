#find patients who may have efficacy to niclosamide

setwd("~/Documents/stanford/hcc/data/")
cancer = "LIHC"
load(paste("sigs", '/', cancer, '_sig_deseq_final.RData', sep='')) #tcga_deseq LIHC_sig_deseq_QC
res$GeneID = sapply((res$id), function(id){unlist(strsplit(id, "\\|"))[2]})
res$symbol = sapply((res$id), function(id){unlist(strsplit(id, "\\|"))[1]})
dz_signature = subset(res, !is.na(padj) & !is.na(id) & id !='?' & padj < 1E-12 & abs(log2FoldChange) > 2 & abs(log2FoldChange) != Inf )
dz_sig = subset(dz_signature, select=c("GeneID", "log2FoldChange"))

data2 = read.csv("microarray//Group3/data2.csv")
data2 = subset(data2, fail.count>=0)
drug_fc = aggregate(Test.Control.fc ~ GENE_ID, data2, mean)
names(drug_fc) = c( "GeneID", "fc")

drug_sample_all = data.frame()

#external matched data
load("~/Documents/stanford/hcc/code/meta_analysis2/GSE54236_GPL6480_matched_expr.RData")
drug_sample = merge(drug_fc, samples_expr, by="GeneID")
drug_sample = subset(drug_sample, GeneID %in% dz_sig$GeneID)
cor(drug_sample, method="spearman")[2,]
p_values = sapply(3:ncol(drug_sample), function(id){
  test = cor.test(drug_sample$fc, drug_sample[,id], method="spearman")
  list(test$p.value, test$estimate)
})
q_values = p.adjust(p_values[1,], "fdr")
length(q_values)
sum(q_values<0.05)

drug_sample_all = rbind(drug_sample_all, data.frame(source = "GSE54236_GPL6480", sample = colnames(drug_sample)[3:ncol(drug_sample)], 
                                                    p = as.vector(unlist(p_values[1,])), q=q_values, rho=as.vector(unlist(p_values[2,]))))

load("~/Documents/stanford/hcc/code/meta_analysis2/GSE14520_GPL3921_matched_expr.RData")
drug_sample = merge(drug_fc, samples_expr, by="GeneID")
drug_sample = subset(drug_sample, GeneID %in% dz_sig$GeneID)
cor(drug_sample, method="spearman")[2,]
p_values = sapply(3:ncol(drug_sample), function(id){
  test = cor.test(drug_sample$fc, drug_sample[,id], method="spearman")
  list(test$p.value, test$estimate)
})
q_values = p.adjust(p_values[1,], "fdr")
length(q_values)
sum(q_values<0.05)

drug_sample_all = rbind(drug_sample_all, data.frame(source = "GSE14520_GPL3921", sample = colnames(drug_sample)[3:ncol(drug_sample)], 
                                                    p = as.vector(unlist(p_values[1,])), q=q_values, rho=as.vector(unlist(p_values[2,]))))



load("~/Documents/stanford/hcc/code/meta_analysis2/GSE14520_GPL571_matched_expr.RData")
drug_sample = merge(drug_fc, samples_expr, by="GeneID")
drug_sample = subset(drug_sample, GeneID %in% dz_sig$GeneID)
cor(drug_sample, method="spearman")[2,]
p_values = sapply(3:ncol(drug_sample), function(id){
  test = cor.test(drug_sample$fc, drug_sample[,id], method="spearman")
  list(test$p.value, test$estimate)
})
q_values = p.adjust(p_values[1,], "fdr")
length(q_values)
sum(q_values<0.05)

drug_sample_all = rbind(drug_sample_all, data.frame(source = "GSE14520_GPL571", sample = colnames(drug_sample)[3:ncol(drug_sample)], 
                                                    p = as.vector(unlist(p_values[1,])), q=q_values, rho=as.vector(unlist(p_values[2,]))))


load("tcga_matched_patient_expr.RData")
drug_sample = merge(drug_fc, patient_sigs, by="GeneID")
drug_sample = subset(drug_sample, GeneID %in% dz_sig$GeneID)
cor(drug_sample, method="spearman")[2,]
p_values = sapply(3:ncol(drug_sample), function(id){
  test = cor.test(drug_sample$fc, drug_sample[,id], method="spearman")
  list(test$p.value, test$estimate)
})
q_values = p.adjust(p_values[1,], "fdr")
length(q_values)
sum(q_values<0.05)

drug_sample_all = rbind(drug_sample_all, data.frame(source = "TCGA", sample = colnames(drug_sample)[3:ncol(drug_sample)], 
                                                    p = as.vector(unlist(p_values[1,])), q=q_values, rho=as.vector(unlist(p_values[2,]))))


drug_sample_all$reversed = "no"
drug_sample_all$reversed[drug_sample_all$q < 0.01 & drug_sample_all$rho<0] = "yes"

library(ggplot2)
pdf("patient_reversed_by_NEN.pdf")
  ggplot(drug_sample_all, aes(source, fill=reversed)) +  theme_bw() +
    geom_bar() + 
    theme(axis.text.x = element_text(size=17, angle = 45, hjust = 1), axis.text.y = element_text(size=17), axis.title.y = element_text(size=18)) +
    xlab("") + 
    ylab("Patient counts")
dev.off()
  
