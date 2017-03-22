#this pipeline is used to refine disease genes that will be used in drug predictions.
#two different sets of disease genes will be used to query CMap and LINCS
#disease gene signatures are essential in the prediction. External sets should be used to validate disease signatures.

cancer = "LIHC"

#choose thresholds that lead to the best separation of tumor/non-tumors using external sets from GEO
#need to reload processed data
source("../code/disease_sig/choose_threshold_for_signatures.R")

#choose thresholds for LINCS 
source("../code/disease_sig/choose_threshold_for_reduced_signatures.R")

#validate disease signatures using the optimal thresholds
source("../code/disease_sig/validate_signatures_external.R")

#validate reduced disease signatures using the optimal thresholds
source("../code/disease_sig/validate_reduced_signatures_external.R")

#validate disease signatures using TCGA 
source("../code/disease_sig/validate_disease_sig.R")

#validate reduced disease signatures using TCGA
source("../code/disease_sig/validate_reduced_disease_sig.R")

