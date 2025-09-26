## Script to compute SAFE-LD starting from summary stats 
# NB data should be provided to this function already cleaned (no NAs) and dimension-matched

## The function takes in input: 
# 1) betas: these should have individuals on the rows and phenotypes on the columns, so it will be of size N x P
# 2) std_err: these should have individuals on the rows and phenotypes on the columns, so it will be of size N x P
# 3) p_values: these should have individuals on the rows and phenotypes on the columns, so it will be of size N x P
# 4) threshold: this is the threshold used to remove SNPs with effect 
# 5) save_RDS: if TRUE, saves SAFE-LD computed from z_scores and from betas in an RDS file called "SAFE-LDss.rds"
# 6) snp_ids: vector containing snp_ids

## The function provides in output: 
# a list containing two objects
# 1) SAFE_LDss computed using z-scores 
# 1) SAFE_LDss computed using betas 

library(data.table)
library(WGCNA)
library(Rfast)

SAFE_LDss <- function (betas, std_errs, p_values, threshold, save_RDS, snp_ids) {
  
  # First we perform standard checks 

  # 1) Check that betas, std_err and p-values all have the same dimension, otherwise
  # there is something wrong with the data input
  if (!identical(dim(betas), dim(std_errs)) || !identical(dim(std_errs), dim(p_values))) {
    stop("Dimensions of betas, std_errs and p_values are not compatible. All three
         should have the same dimension and they should be formatted so that individuals
         are on the rows and phenotypes on the columns. ")
  } 
  
  # 2) Now we check if NAs are present 
  if (sum(is.na(betas)) != 0) {
    stop("Betas provided contain NAs")
  } else if (sum(is.na(std_errs)) != 0) {
    stop("Standard errors provided contain NAs")
  } else if (sum(is.na(p_values)) != 0) {
    stop("P-values provided contain NAs")
  }
  
  # 3) Give warning if number of individuals and number of phenotypes are very different
  # as this may impact quality of fine-mapping
  n = dim(betas)[1]
  p = dim(betas)[2]
  
  if (p < n/2) {
    message("WARNING: Number of phenotypes is less than half the number of individuals 
            (P < N/2). Care should be taken as this could lead to low power
            and coverage in the application of fine-mapping.")
  }
  
  # 4) Give error if length of snp_ids does not correspond to number of individuals 
  if (n != length(snp_ids)) {
    stop("Length of SNP IDs provided does not match number of individuals.")
  }
  
  # Now we start with SAFE-LDss
  
  # The first step is to perform thresholding based on the threshold provided by the user 
  # We select only those which do not have an effect
  idxs_to_keep <- which(colMins(p_values) > threshold)
  betas_reduced <- betas[,idxs_to_keep]
  #num_traits_kept <- ncol(data.frame(betas_reduced))
  #percentage_phenos <- (num_traits_kept/pheno_sizes[j]) * 100
  se_reduced <- std_errs[,idxs_to_keep]
  
  zscores <- matrix(0, nrow = dim(betas_reduced)[1], dim(betas_reduced)[2])
  
  for (i in 1:ncol(betas_reduced)) {
    zscores[,i] <- betas_reduced[,i] / se_reduced[,i]
  }
  
  # Now we compute SAFE-LD
  betas_cor <- WGCNA::cor(t(betas_reduced))
  rownames(betas_cor) <- snp_ids
  colnames(betas_cor) <- snp_ids
  zscore_cor <- WGCNA::cor(t(zscores)) 
  rownames(zscore_cor) <- snp_ids
  colnames(zscore_cor) <- snp_ids
  
  # Put desired outputs in a list
  output <- list(
    `zscore_cor` = zscore_cor,
    `betas_cor` = betas_cor
  )
  
  # Save as an RDS
  if (save_RDS == TRUE) {
    saveRDS(output, "SAFE-LDss.rds")
  }
  
  return(output)
  
}


                  
