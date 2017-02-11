# run cmap in lincs

args <- commandArgs(trailingOnly=T)
subset_comparison_id <-  args[1]
analysis_id <- args[2]
dz_sig_path <- args[3]
landmark <- args[4] #1 means using landmark gene, 0 means using all genes

if (length(args)<4){
	print ("please input 3 arguments (subset comparision id, analysis id, disease signature path, landmark)")
	q()
}

cmap_score_old <- function(sig_up, sig_down, drug_signature) {
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

cmap_score <- function(sig_up, sig_down, drug_signature) {
	
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

	if(num_tags_up > 1 && num_tags_down > 1) {
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


if (landmark == 1){
  landmark_data = read.csv("raw/lincs/lincs_landmark.csv")
  gene.list = landmark_data$gene_id
  load("raw/lincs/lincs_signatures_cmpd_landmark.RData")
  
}else{
  landmark_data = read.csv("raw/lincs/probe_id_info.csv")
  gene.list = landmark_data$gene_id
  load("raw/lincs/lincs_signatures_cmpd_hcc_all.RData")
}

lincs_sig_info = read.csv("raw/lincs/lincs_sig_info.csv")
lincs_sig_info = subset(lincs_sig_info, id %in% colnames(lincs_signatures))
#remove duplicate instances
lincs_sig_info = lincs_sig_info[!duplicated(lincs_sig_info$id),]

sig.ids = lincs_sig_info$id

# Load in the signature for the subset_comparison_id
dz_signature <- read.table(dz_sig_path,header=T,sep="\t")
dz_signature <- subset(dz_signature, GeneID %in% gene.list)

dz_genes_up <- subset(dz_signature,up_down=="up",select="GeneID")
dz_genes_down <- subset(dz_signature,up_down=="down",select="GeneID")

max_gene_size = 150
#only select 100 genes
if (nrow(dz_genes_up)> max_gene_size){
	dz_genes_up <- data.frame(GeneID=dz_genes_up[1:max_gene_size,])
}
if (nrow(dz_genes_down)> max_gene_size){
	dz_genes_down <- data.frame(GeneID=dz_genes_down[1:max_gene_size,])
}

print ((nrow(dz_genes_up)+nrow(dz_genes_down)))

N_PERMUTATIONS <- 10000 #default 100000
random_sig_ids = sample(colnames(lincs_signatures),N_PERMUTATIONS,replace=T)
count = 0
random_cmap_scores = NULL
for (exp_id in random_sig_ids){
  count = count + 1
  print(count)
  if (landmark ==1){
    cmap_exp_signature <- data.frame(gene.list,  rank(-1 * lincs_signatures[, as.character(exp_id)], ties.method="random"))    
  }else{
    cmap_exp_signature <- data.frame(gene.list,  lincs_signatures[, as.character(exp_id)])    
  }
  colnames(cmap_exp_signature) <- c("ids","rank")
  
  random_input_signature_genes <- sample(gene.list, (nrow(dz_genes_up)+nrow(dz_genes_down)))
  rand_dz_gene_up <- data.frame(GeneID=random_input_signature_genes[1:nrow(dz_genes_up)])
  rand_dz_gene_down <- data.frame(GeneID=random_input_signature_genes[(nrow(dz_genes_up)+1):length(random_input_signature_genes)])
  random_cmap_scores = c(random_cmap_scores, cmap_score(rand_dz_gene_up,rand_dz_gene_down,cmap_exp_signature))
}

rand_cmap_scores = random_cmap_scores 
save(rand_cmap_scores,file=paste(subset_comparison_id, '/drug/',subset_comparison_id,'_',analysis_id,"_", landmark, '_randoms.RData',sep=""))

