# Computes the ConnectivityMap score based on Marina's implementation of the KS-based CMap algorithm
# $1 = subset_comparison_id, $2 = analysis_id
library(qvalue)

source("/projects/reposition/code/compute/new_set/db_conn.R")
args <- commandArgs(trailingOnly=T)
subset_comparison_id <- args[1]
analysis_id <-  args[2]
dz_sig_path <- args[3]

if (length(args)<3){
        print ("please input 3 arguments (subset comparision id, analysis id, disease signature path)")
        q()
}

cmap_score <- function(sig_up, sig_down, drug_signature) {
	# Note. I think that creating the anonymous functions in each iteration of the sapply's below is slowing things down. Predefine them eventually.
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

	# this case, no down regulated gene
	if(num_tags_up > 1 || num_tags_down > 1) {
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
		
		if (sum(sign(c(ks_down,ks_up))) == 0) {
			connectivity_score <- ks_up - ks_down # different signs
		}
	}
	return(connectivity_score)
}

# cmap_score <- function(sig_up, sig_down, drug_signature) {
# 	
# 	num_genes = nrow(drug_signature)
# 
# 	ks_up <- 0
# 	ks_down <- 0
# 
# 	connectivity_score <- 0
# 
# 	drug_signature[,"rank"] <- rank(drug_signature[,"rank"])
# 
# 	# Merge the drug signature with the disease signature by GeneID
# 	temp_up <- merge(drug_signature, sig_up, by.x = "ids", by.y = 1)
# 	temp_down <- merge(drug_signature, sig_down, by.x = "ids", by.y = 1)
# 
# 	temp_up <- sort(temp_up$rank)
# 	temp_down <- sort(temp_down$rank)
# 
# 	num_tags_up <- length(temp_up)
# 	num_tags_down <- length(temp_down)
# 
# 	if(num_tags_up > 1 && num_tags_down > 1) {
# 		
# 		a_up <- 0
# 		b_up <- 0
# 		
# 		for(j in 1:num_tags_up) {
# 			a_iter <- j/num_tags_up - temp_up[j]/num_genes
# 			b_iter <- temp_up[j]/num_genes - (j-1)/num_tags_up
# 
# 			if(a_iter > a_up) {
# 				a_up <- a_iter
# 			}
# 			if(b_iter > b_up) {
# 				b_up <- b_iter
# 			}
# 		}
# 
# 		if(a_up > b_up) {
# 			ks_up <- a_up
# 		} else {
# 			ks_up <- -b_up
# 		}
# 
# 		a_down <- 0
# 		b_down <- 0
# 
# 		for(j in 1:num_tags_down) {
# 			a_iter <- j/num_tags_down - temp_down[j]/num_genes
# 			b_iter <- temp_down[j]/num_genes - (j-1)/num_tags_down
# 
# 			if(a_iter > a_down) {
# 				a_down <- a_iter
# 			}
# 
# 			if(b_iter > b_down) {
# 				b_down <- b_iter
# 			}
# 		}
# 
# 		if(a_down > b_down) {
# 			ks_down <- a_down
# 		} else {
# 			ks_down <- -b_down
# 		}
# 		
# 		if((ks_down > 0 && ks_up > 0) || (ks_down < 0 && ks_up < 0)) {
# 			connectivity_score <- 0
# 		} else {
# 			connectivity_score <- ks_up - ks_down
# 		}
# 	}
# 	return(connectivity_score)
# }

# So the strategy will be. filter/sort by fdr, filter/sort by p-value, sort by fold-change up to 100 genes
# If we have an experiment where the fdr does not meet our cutoff, just use p-value and fold-change cutoff to make a set of 10 genes
get_dz_signature <- function(subset_comparison_id,analysis_id) {
	p_thresh <- 0.01
	q_thresh <- 0.05
	res <- dbSendQuery(
		db_con,
		paste(
			"SELECT * FROM analysis_results WHERE analysis_id = ",analysis_id," AND subset_comparison_id = ",subset_comparison_id," AND GeneID IS NOT NULL",
			sep=""
		)
	)
	dz_siggenes <- fetch(res,n=-1)
	
	# dz_up_genes <- subset(dz_siggenes,comment=="DZUP" & sig_value <= p_thresh,select=c('GeneID','delta_value','sig_value','fdr_value'))
	# dz_down_genes <- subset(dz_siggenes,comment=="DZDOWN" & sig_value <= p_thresh,select=c('GeneID','delta_value','sig_value','fdr_value'))
	# Changing to SAMR to compare with marina
	dz_up_genes <- subset(dz_siggenes,comment=="DZUP" & fdr_value <= q_thresh,select=c('GeneID','delta_value','sig_value','fdr_value'))
	dz_down_genes <- subset(dz_siggenes,comment=="DZDOWN" & fdr_value <= q_thresh,select=c('GeneID','delta_value','sig_value','fdr_value'))
	
	
	dz_up_sorted <- dz_up_genes[order(dz_up_genes$fdr_value,-dz_up_genes$delta_value),]
	dz_down_sorted <- dz_down_genes[order(dz_down_genes$fdr_value,dz_down_genes$delta_value),]
	
	return(
		list(
			up.sig = dz_up_sorted[1:(min(nrow(dz_up_sorted),100)),],
			down.sig = dz_down_sorted[1:(min(nrow(dz_down_sorted),100)),]
		)
	)
}
#load gist genes
dz_signature <- read.table(dz_sig_path,header=T,sep="\t")
dz_genes_up <- subset(dz_signature,up_down=="up",select="GeneID")
dz_genes_down <- subset(dz_signature,up_down=="down",select="GeneID")

#only choose the top 100 genes
if (nrow(dz_genes_up)>100){
#        dz_genes_up <- data.frame(GeneID= dz_genes_up[1:100,])
}
if (nrow(dz_genes_down)>100){
#        dz_genes_down <- data.frame(GeneID=dz_genes_down[1:100,])
}

# Load in the signature for the subset_comparison_id
#dz_signature <- get_dz_signature(subset_comparison_id,analysis_id)
#dz_genes_up <- subset(dz_signature$up.sig,select=GeneID)
#dz_genes_down <- subset(dz_signature$down.sig,select=GeneID)

#load('/home/jdudley/scratch/geneid_processed_data_all.RData')
load('/projects/reposition/data/cmap_signatures_updated.RData')
gene_list <- subset(cmap_signatures,select=1)
cmap_signatures <- cmap_signatures[,2:ncol(cmap_signatures)] # Shouldn't overwrite, but worried about memory usage
dz_cmap_scores <- sapply(1:ncol(cmap_signatures),function(exp_id) {
	print(paste("Computing score for disease against cmap_experiment_id =",exp_id))
	cmap_exp_signature <- cbind(gene_list,subset(cmap_signatures,select=exp_id))
	colnames(cmap_exp_signature) <- c("ids","rank")
	cmap_score(dz_genes_up,dz_genes_down,cmap_exp_signature)
})

# Compute the significance
print("Loading random scores for subset")
load(paste('/projects/reposition/cmap_randoms/',subset_comparison_id,'_',analysis_id,'_randoms.RData',sep=""))
random_scores <- unlist(rand_cmap_scores)
# Frequency-based p-value using absolute scores from sampling distribution to approximate two-tailed p-value
print("COMPUTING p-values")
p_values <- sapply(dz_cmap_scores,function(score) {
	length(which(abs(random_scores) >= abs(score))) / length(random_scores)
})
print("COMPUTING q-values")
q_values <- qvalue(p_values)$qvalues

# Save the results to the database
#dbSendQuery(db_con,"BEGIN")

	sapply(1:length(dz_cmap_scores),function(exp_id) {
		dbSendQuery(db_con,paste("INSERT INTO proj_repositioning.cmap_predictions_temp (experiment_id,subset_comparison_id,cmap_score,p_value,q_value,analysis_id) VALUES (",exp_id,",","\"",subset_comparison_id,"\"",",",dz_cmap_scores[exp_id],",",p_values[exp_id],",",q_values[exp_id],",",analysis_id,")",sep=""))
	})

#dbSendQuery(db_con,"COMMIT")
