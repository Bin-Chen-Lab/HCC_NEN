#need to set up connection to database, db_conn.R

#predict drugs using CMap database
Rscript predict_drugs_cmap.R subset_comparison_id analysis_id signature_path

#predict drugs using lincs database
Rscript predict_drugs_lincs.R subset_comparison_id analysis_id signature_path flag
flag = 1: using landmark genes
flag = 0: using all genes 

