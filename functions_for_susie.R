library(genio)
library(data.table)
library(susieR)
library(pgenlibr)
library(WGCNA)
library(stringr)
library(dplyr)

source("/home/giulia.desanctis/ld_sim_git/ld-simulation-giulia/bin/funs_locus_breaker_cojo_finemap_all_at_once.R")


susie_ht <- function(beta, se, R, var_y = 1, L = 10, N, estimate_residual_variance, max_iter){
  
  fitted_rss <- susieR::susie_rss(
    bhat = beta, 
    shat = se, 
    n = N, 
    R = R, 
    var_y = 1, 
    L = L,
    estimate_residual_variance = estimate_residual_variance,
    max_iter = max_iter
    
  )
  return(fitted_rss)
  
}

do_susie <- function(sim, R_to_use, est_res_var, name_of_LD, model_no, num_ind, max_iter=2000, num_phenos=NULL) {
  susie.cs=list()
  overlap_snps <- intersect(sim$sum_stats$ID,colnames(R_to_use))
  sim$sum_stats <- sim$sum_stats %>% filter(ID %in% overlap_snps)
  R_to_use <- R_to_use[overlap_snps,overlap_snps]
  
  susie.cs <- vector("list", length(5))
  
  for(i in 1:length(sim$effect_list)){
    
    dataset=dataset.munge_hor( sumstats.file = sim$sum_stats
                               ,chr.lab ="#CHROM"
                               ,snp.lab = "ID"
                               ,pos.lab = "POS"
                               ,a1.lab = "A1"
                               ,a0.lab = "REF"
                               ,beta.lab = paste0("P",i,"_BETA")
                               ,se.lab = paste0("P",i,"_SE")
                               ,freq.lab = "A1_FREQ"
                               ,sdY = 1
                               ,type = "quant"
                               ,n.lab = "N")
    chr = unique(dataset$CHR)
    dataset$phenotype_id="Pippo"
    all.dataset=dataset.align(dataset = dataset)
    all.dataset$SNP=all.dataset$snp_original
    
    susie.res <- try(susie_ht(
      beta = sim$sum_stats[[paste0("P",i,"_BETA")]], # sim$s_stat$P1_BETA
      se = sim$sum_stats[[paste0("P",i,"_SE")]], # sim$s_stat$P1_BETA
      R = R_to_use,
      var_y = 1,
      L = 10,
      N = num_ind,
      estimate_residual_variance = est_res_var, 
      max_iter = max_iter
    ), outFile = getOption("try.outFile", default = stderr()))
    
    if (inherits(susie.res, "try-error")) {
      warning(paste("susie_ht failed:", as.character(susie.res)))
    }
    
    susie.cs[[i]] = susie.res
    print(sprintf("iteration: %f", i))
    
  }
  
  if(name_of_LD == "simulated") {
  saveRDS(
    susie.cs,
    file=sprintf("susie_model_%d_with_%s_LD_%d_phenos_%d_individuals.rds", model_no, as.character(name_of_LD), num_phenos, num_ind)
  )
  } else {
    saveRDS(
      susie.cs,
      file=sprintf("susie_model_%d_with_%s_LD.rds", model_no, as.character(name_of_LD))
    )
  }
}

compute_coverage <- function(effect_list, susie_out) {
  coverage_num <- integer(length(susie_out))
  coverage_den <- integer(length(susie_out))
  
  for (i in 1:length(susie_out)) {
    susie <- susie_out[[i]]
    
    if (inherits(susie, "try-error")) {
      warning(sprintf("Skipping i = %d due to susie error", i))
      coverage_num[i] <- 0
      coverage_den[i] <- 0
      next
    }
    
    all_snps <- colnames(susie$alpha)
    
    cs_sets <- susie$sets$cs
    if (is.null(cs_sets)) {
      coverage_num[i] <- 0
      coverage_den[i] <- 0
      next
    }
    
    true_causal_snps <- effect_list[[i]]$snps
    coverage_den[i] <- length(cs_sets)
    
    overlap_counts <- sapply(cs_sets, function(cs) {
      any(all_snps[unlist(cs)] %in% true_causal_snps)
    })
    
    coverage_num[i] <- sum(overlap_counts)
  }
  
  return(list(
    coverage_num = coverage_num,
    coverage_den = coverage_den
  ))
}

compute_power <- function(effect_list, susie_out) {
  power_num <- integer(length(susie_out))
  power_den <- integer(length(susie_out))
  
  for (i in seq_along(susie_out)) {
    
    if (inherits(susie, "try-error")) {
      warning(sprintf("Skipping i = %d due to susie error", i))
      coverage_num[i] <- 0
      coverage_den[i] <- 0
      next
    }
    
    susie <- susie_out[[i]]
    
    cs_sets <- susie$sets$cs
    if (is.null(cs_sets)) {
      power_num[i] <- 0
      power_den[i] <- length(effect_list)
      next
    }
    
    snp_names <- colnames(susie$mu)
    power_den[i] <- length(effect_list[[i]])
    
    snps_in_cs <- unique(unlist(lapply(cs_sets, function(cs) snp_names[unlist(cs)])))
    
    power_num[i] <- sum(effect_list[[i]]$snps %in% snps_in_cs)
  }
  
  return(list(
    power_num = power_num,
    power_den = power_den
  ))
}

