#this pipeline is used to predict drug hits 

cancer = "LIHC"

#prepare signatures for cmap and lincs
#source("../code/drug_prediction/build_dz_signature.R")

#predict cmap drugs. Take a couple of hours
cmd = paste("Rscript ../code/drug_prediction/predict_drugs_cmap.R", cancer, "cmap", paste(cancer, "/drug/dz_signature_cmap.txt", sep=""))
system(cmd)

#predict lincs drugs. Take a few of hours
#cmd = paste("Rscript ../code/drug_prediction/predict_drugs_lincs.R", cancer, "lincs", paste(cancer, "/drug/dz_signature_lincs.txt", sep=""), 1)

cmd = paste("Rscript ../code/drug_prediction/predict_drugs_lincs.R", cancer, "lincs", paste(cancer, "/drug/dz_signature_cmap.txt", sep=""), 0)

system(cmd)

#analyze predictions
source("../code/drug_prediction/drug_repurpose.R")