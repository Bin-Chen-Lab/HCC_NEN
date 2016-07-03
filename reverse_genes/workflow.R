#this pipeline is used to examine the correlation between disease signatures and drug signatures derived from in vitro/in vivo,
#as well as to find disease genes reversed by niclosamide/NEN

cancer = "LIHC"

#compute correlation and identify reverse genes
source("../code/reverse_genes/reverse_genes.R")

#