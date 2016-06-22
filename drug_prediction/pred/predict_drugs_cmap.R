#predict using cmap data
#the CMAP score computation is based on the Marina's code (PMC3502016)

args <- commandArgs(trailingOnly=T)

subset_comparison_id <- args[1]
analysis_id <-  args[2]
dz_sig_path <- args[3]

if (length(args)<3){
  print ("please input 3 arguments (subset comparision id, analysis id, disease signature path)")
  q()
}

#compute random scores, used for signficance analysis
cmd = paste("Rscript compute_cmap_randoms_new_set.R", subset_comparison_id, analysis_id, dz_sig_path)
#system(cmd)

#compute cmap scores 
cmd = paste("Rscript compute_connectivity_score_new_set.R", subset_comparison_id, analysis_id, dz_sig_path)
system(cmd)
