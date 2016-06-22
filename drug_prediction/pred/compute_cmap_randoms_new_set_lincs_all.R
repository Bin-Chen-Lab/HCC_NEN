# run cmap in lincs

#source("/projects/reposition/code/compute/new_set/db_conn.R")

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


get.cmpd.sig.file <- function(id,  landmark=F){
  #reading sig from file is much faster
  if (landmark == T){
    rows = 978
  }else{
    rows = -1
  }
  sigs = read.table(paste("/home/binchen1/scratch/lincs/data/sig/", id, ".txt", sep=""), nrows=rows)
  return(sigs[,1])
}

get.sig.ids <- function(){
  cmpd_sigs = read.delim("/home/binchen1/proj/cmap/data/cmpd_sig_info_gold.txt", header=T, sep="\t", stringsAsFactors=F)
  return(unique(cmpd_sigs$id))
}

get.sigs <- function(id){
  #rank, the higher zcore was ranked in the bottom, but cmap score assumes it is in the top, so reverse the rank here
  return(-1 * scan(paste("/home/binchen1/scratch/lincs/data/sig_rank/", id, ".txt", sep="")))
}

get.sigs.landmark <- function(id){
  # return(1 * scan(paste("/home/binchen1/scratch/lincs/data/sig_rank/", id, ".txt", sep="")))
  
  sig.raw = (get.cmpd.sig.file(id, T))
  
  return(rank(-sig.raw, ties.method="random"))
}

get.lincs.gene.list <- function(con){
  lincs_gene <- dbReadTable(con, "probe_id_info")
  
  return(lincs_gene$gene_id[1:978])
}


build.connection <- function(){
  library("RMySQL")
  mysql_drvr <-dbDriver("MySQL")
  con <- dbConnect(mysql_drvr, group="client",host="buttelab-db1.stanford.edu",dbname="proj_lincs")
}

get.gene.list <- function(con){
  return(scan(paste("/home/binchen1/scratch/lincs/data/sig_rank/gene_id.txt", sep="")))
}

con = build.connection()
if (landmark == 1){
  gene.list = get.lincs.gene.list(con)
}else{
  gene.list = get.gene.list()  
}

sig.ids = get.sig.ids()

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

#load('/home/binchen1/scratch/geneid_processed_data_all.RData')
#load('/projects/reposition/data/geneid_processed_data_all.RData')
#gene_list <- subset(cmap_signatures,select=1)
#cmap_signatures <- cmap_signatures[,2:ncol(cmap_signatures)] # Shouldn't overwrite, but worried about memory usage
N_PERMUTATIONS <- 100000 #default 100000


# So Basically, given a signature with this num_tags_up and num_tags_down, what is the null cmap_score distribution?
rand_cmap_scores <- sapply(sample(sig.ids,N_PERMUTATIONS,replace=T),function(exp_id) {
	print(paste("Computing score for disease against cmap_experiment_id =",exp_id))
	if (landmark == 1){
	  cmap_exp_signature <- cbind(gene.list,  get.sigs.landmark(exp_id))    
	}else{
	  cmap_exp_signature <- cbind(gene.list,  get.sigs(exp_id))    
	}
	colnames(cmap_exp_signature) <- c("ids","rank")
	random_input_signature_genes <- sample(gene.list, (nrow(dz_genes_up)+nrow(dz_genes_down)))
	rand_dz_gene_up <- data.frame(GeneID=random_input_signature_genes[1:nrow(dz_genes_up)])
	rand_dz_gene_down <- data.frame(GeneID=random_input_signature_genes[(nrow(dz_genes_up)+1):length(random_input_signature_genes)])
	cmap_score(rand_dz_gene_up,rand_dz_gene_down,cmap_exp_signature)
},simplify=F)

# Relational DB is too painful for stuff like this. Saving to disk. Need a hash DB.
save(rand_cmap_scores,file=paste('/projects/reposition/cmap_randoms/',subset_comparison_id,'_',analysis_id,'_', landmark, '_randoms.RData',sep=""))

# sapply(1:length(rand_cmap_scores),function(exp_id) {
# 	dbSendQuery(db_con,paste("INSERT INTO proj_repositioning.cmap_randoms (experiment_id,rand_round,subset_comparison_id,cmap_score) VALUES (",exp_id,",",N_PERMUTATIONS,",",subset_comparison_id,",",rand_cmap_scores[exp_id],")",sep=""))
# })
