## Script to compute SAFE-LD starting from summary stats 

#1) are there infinite betas
#2) infinite SE
#3) remove missing data 
#import gwas, do basic checks (no. of rows, is there a beta and se or z-score, no infinite )
#munch sum stats —> harmonizes gwas data that cleans up gwas 
#have code print warnings  if P and N are unblanced  post-filtering 


## The function takes in input: 
# 1) betas: these should have individuals on the rows and phenotypes on the columns 
# 2) std_err: these should have individuals on the rows and phenotypes on the columns 
# 3) p_values: these should have individuals on the rows and phenotypes on the columns
# 4) threshold: this is the threshold used to remove SNPs with effect 
# 5) save_RDS: if TRUE, saves SAFE-LD with z_scores and with betas in an RDS file 
# 6) snp_ids: vector containing snp_ids

## The function provides in output: 
## a list containing SAFE_LDss using betas and SAFE_LDss using z-scores

library(data.table)
library(WGCNA)
library(Rfast)
library(MungeSumstats)

SAFE_LDss <- function (betas, std_errs, p_values, threshold, save_RDS, snp_ids) {
  
  # First we check that betas, std_err and p-values all have the same dimension, otherwise
  # there is something wrong with the data input

  # Compare the dimensions
  if (!identical(dim(betas), dim(std_errs)) || !identical(dim(std_errs), dim(p_values))) {
    stop("Dimensions of betas, std_errs and p_values are not compatible. All three
         should have the same dimension and they should be formatted so that individuals
         are on the rows and phenotypes on the columns. ")
  } 
  
  # Now we check if NAs are present 
  if (sum(is.na(betas)) != 0) {
    stop("Betas provided contain NAs")
  } else if (sum(is.na(std_errs)) != 0) {
    stop("Standard errors provided contain NAs")
  } else if (sum(is.na(p_values)) != 0) {
    stop("P-values provided contain NAs")
  }
  
  # Give warning if number of individuals and number of phenotypes are very different
  n = dim(betas)[1]
  p = dim(betas)[2]
  
  if (p < n/2) {
    message("WARNING: Number of phenotypes is less than half the number of individuals 
            (P < N/2). Care should be taken as this could lead to low power
            and coverage in the application of fine-mapping.")
  }
  
  # Now we perform thresholding 
  # Selecting only those which do not have an effect
  idxs_to_keep <- which(colMins(p_values) > threshold)
  betas_reduced <- betas[,idxs_to_keep]
  num_traits_kept <- ncol(data.frame(betas_reduced))
  #percentage_phenos <- (num_traits_kept/pheno_sizes[j]) * 100
  se_reduced <- std_errs[,idxs_to_keep]
  
  zscores <- matrix(0, nrow = dim(betas_reduced)[1], dim(betas_reduced)[2])
  zscores <- as.data.frame(zscores)
  
  for (i in 1:ncol(betas_reduced)) {
    zscores[[i]] <- betas_reduced[[i]] / se_reduced[[i]]
  }
  
  # Now we compute SAFE-LD
  
  betas_cor <- WGCNA::cor(t(betas_reduced))
  #rownames(betas_cor) <- snp_ids
  #colnames(betas_cor) <- snp_ids
  zscore_cor <- WGCNA::cor(t(zscores)) 
  #rownames(zscore_cor) <- snp_ids
  #colnames(zscore_cor) <- snp_ids
  
  
  output <- list(
    `betas_cor` = betas_cor,
    `zscore_cor` = zscore_cor
  )
  
  if (save_RDS == TRUE) {
    saveRDS(output, "SAFE-LDss")
  }
  
  return(output)
  
}


                  
