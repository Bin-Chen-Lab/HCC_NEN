#this pipeline can be used for any other cancers studied in TCGA.
#author: Bin Chen (July 2014)

###################
#parameters
##################
cancer =  'LIHC' 
cancer_name = 'Hepatocellular Carcinoma' #used for label

#T: RNA-Seq data from GDAC; #F: manually downloaded from TCGA website
data_from_gdac = T 

#varying5k: 5000 varying genes; also support other gene sets (sigs: only DE genes from the meta-analysis;  metabolism: metabolism realted genes
# meta_sigs: DE genes from meta analysis, #varying5k_tumor)
comparison_gene_set = "varying5k"
#number of varying genes used to correlate tumors and cell lines
num_varying_genes = 5000

#########
#MAIN code
#########

#download tumors from GDAC
source("../code/tumor_cell_line/download_from_gdac.R")

#download clincal data and mrna meta data from TCGA
#####need to go to TCGA website to download all the required data; 

#compute random tumors and cell line correlations. Take a few minutes. 
cutoff = 9.34E-06 # #comment out this line if you want to use the random samples to correct the p values. 
#source("../code/tumor_cell_line/correct_by_expo.R")

#compute tumor vs cell line correlations
source("../code/tumor_cell_line/compute_tumor_cell_line_cor_update.R")

#analyze correlations and select cell lines
source("../code/tumor_cell_line/select_cell_lines_update.R")

#compute differentially expressed genes between tumors and non-tumors
source("../code/tumor_cell_line/compute_disease_signatures.R")

#compute differentially expressed genes between tumors and cell lines
#source("../code/tumor_cell_line/comput_tumor_cell_line_diff.R")

#assess correlation between tumors from external datasets and cell lines 
source("../code/tumor_cell_line/external_dataset_cor_stats.R")
