LIHC_sig_sample_cutoff = 9.34E-06 #from expo geo correction

#compute ratio
gpl = "GPL8177"
data = read.csv(paste("tumor_cell_line/tumor_cell_allLIHC", gpl, ".csv", sep=""))

sample_class = read.table("../code/meta_analysis2/meta_input.txt", sep="\t", header=T)
sample_class_tumor = subset(sample_class, ClassCode ==2, select = c("GSM"))
sample_cor = aggregate(outlier ~ sample_id, data, median)
sample_cor_tumor = subset(sample_cor, sample_id %in% sample_class_tumor$GSM)

dim(sample_cor_tumor)
sum(sample_cor_tumor$outlier < LIHC_sig_sample_cutoff)

#compute median; range
gpls = c("GPL6480", "GPL14550", "GPL10558", "GPL10687", "GPL3921", "GPL571")
sample_class = read.table("../code/meta_analysis2/meta_input.txt", sep="\t", header=T)
sample_class_tumor = subset(sample_class, ClassCode ==2, select = c("GSM"))

data = data.frame()
for (gpl in gpls){
  data_gpl = read.csv(paste("tumor_cell_line/tumor_cell_allLIHC", gpl, ".csv", sep=""))
  data = rbind(data, data_gpl)
}
data = subset(data, sample_id %in% sample_class_tumor$GSM & cell_line_name %in% c("PLC/PRF/5", "HuH-7", "Hep G2"))
summary(data$cor)
cell_line_cor = aggregate(cor ~ cell_line_name, data, median)

