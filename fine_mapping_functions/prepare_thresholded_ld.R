args <- commandArgs(trailingOnly = TRUE)
model_number <- as.integer(args[1])

library(genio)
library(data.table)
library(susieR)
library(pgenlibr)
library(WGCNA)
library(stringr)
library(pracma)
library(abind)
library(matrixStats)
library(readr)

beta.calc=function(r2,maf) sqrt(r2/(2*maf*(1-maf)))

prepare_thresholded_ld <- function(
     model.no=1
    ,null.phenotype.path="/project/statgen/models_for_final_simulations/model_%d/20k_pheno_file_model_%d.rds"
    ,effect.phenotype.path = "/scratch/giulia.desanctis/susie_results_method_2/phenotypes/model_%d/20k_pheno_file_model_%d_effect_%f.rds"
    ,pheno_sizes=c(1000)
    ,bed.internal="/ssu/bsssu/giulia_ld_simulations/2025_05_02_brit_iris_genotype_data/ukbb_all_chrs_grch38_maf0.01_30000_random_unrelated_white_british_alpha_sort_alleles_aligned_nodups_for_sim" 
    ,plink2.binary="/ssu/gassu/software/plink/2.00_20211217/plink2"
    ,plink.binary="/ssu/gassu/software/plink/1.90_20210606/plink"
    ,maf.threshold=0.01
    ,chr=6
    ,start.pos=134597	
    ,end.pos=135597			
    ,keep.snps=FALSE
    ,snp.list=NULL
    ,subsample_individuals=TRUE
    ,individuals="~/individual_lists/individual_idxs_internal_10k.txt"
    ,num_individuals=1000
    ,b.file.extension="all_indiv"
    ,b.file.extension.subsample="subsample_indiv"
    ,keep.bfile=FALSE
    ,effect.sizes=c(0.05)
    ,dec.phenos.with.effect=0.1
    ,threshold=1e-6
) 
{
  
  print(sprintf("doing model: %d", model.no))
  
  perc_phenos_with_effect = dec.phenos.with.effect * 100
  
  clean_perc <- ifelse(perc_phenos_with_effect %% 1 == 0,
                       sprintf("%d", perc_phenos_with_effect),
                       sprintf("%.1f", perc_phenos_with_effect))
  
  setwd(sprintf("/scratch/giulia.desanctis/susie_results_method_2/%d_individuals/model_%d/%s_perc_phenos_with_effect/effect_%.2f",
          num_individuals, model.no, clean_perc, effect.sizes))
  
  # 1) Using internal bed files 
  
  all_phenos_nulls <- readRDS(sprintf(null.phenotype.path, model.no, model.no))
  all_phenos_effect <- readRDS(sprintf(effect.phenotype.path, model.no, model.no, effect.sizes))
  
  # Now we need to subset the phenotypes
  phenos_of_different_sizes <- list()
  
  for (i in 1:length(pheno_sizes)) {
    num_phenos_with_effect <- dec.phenos.with.effect * pheno_sizes[i]
    num_phenos_nulls <- pheno_sizes[i] - num_phenos_with_effect
    idxs_effect <- sample(3:5002, num_phenos_with_effect, replace = FALSE)
    idxs_nulls <- sample(3:20002, num_phenos_nulls, replace = FALSE)
    phenos_effect <- all_phenos_effect$phenos[,c(idxs_effect)]
    phenos_nulls <- all_phenos_nulls$phenos[,c(1,2,idxs_nulls)]
    phenos <- cbind(phenos_nulls, phenos_effect)
    new_names <- paste0("P", seq_len(pheno_sizes[i]))
    colnames(phenos) <- c("FID", "IID", new_names)
    if (subsample_individuals == TRUE) {
      inds <- read.table(individuals)
      inds_id <- unlist(inds$V1)
      phenos <- phenos[phenos$FID %in% inds_id,]
      rownames(phenos) <- 1:num_individuals
    }
    phenos_of_different_sizes[[i]] <- phenos
  }
  
  # Now we need to calculate the different LD matrices 
  all_summaries = list()
  matrix_names = list()
  r_sim_betas = list()
  r_sim_zscore = list()
  perc_phenos_remaining_post_thresholding = list()
  
  if (subsample_individuals == TRUE) {
    commando2=paste(plink2.binary,"--bfile",bed.internal, "--keep",individuals, "--keep-allele-order","--chr",chr,"--from-kb",start.pos,"--to-kb",end.pos,"--maf",maf.threshold, "--make-bed", "--out", b.file.extension.subsample)
    system(commando2)
    genos = read_plink(b.file.extension.subsample)
  } else {
    commando1=paste(plink2.binary,"--bfile",bed.internal,"--keep-allele-order","--chr",chr,"--from-kb",start.pos,"--to-kb",end.pos,"--maf",maf.threshold,"--make-bed", "--out", b.file.extension)
    system(commando1)
    genos = read_plink(b.file.extension)
  }
  
  for (j in seq(length(phenos_of_different_sizes))) {
    sprintf("number of phenos: %d", pheno_sizes[j])
    fwrite(phenos_of_different_sizes[[j]], file = sprintf("%d_phenotypes.phen", pheno_sizes[j]), row.names = FALSE, col.names = TRUE, quote = FALSE, sep = "\t")
    
    if (subsample_individuals == TRUE) {
      command1=paste(plink2.binary,"--bfile",b.file.extension.subsample, "--pheno", sprintf("%d_phenotypes.phen", pheno_sizes[j])," --glm cols=+a1freq allow-no-covars omit-ref --out",b.file.extension.subsample) 
    } else {
      command1=paste(plink2.binary,"--bfile",b.file.extension, "--pheno", sprintf("%d_phenotypes.phen", pheno_sizes[j])," --glm cols=+a1freq allow-no-covars omit-ref --out",b.file.extension) 
    }
    system(command1)
    
    for (i in colnames(phenos_of_different_sizes[[j]])[-c(1:2)]) {
      if(subsample_individuals == TRUE) {
        file.name = paste0(b.file.extension.subsample, ".", i, ".glm.linear")
      } else {
        file.name = paste0(b.file.extension, ".", i, ".glm.linear")
      }
      
      if(i == colnames(phenos_of_different_sizes[[j]])[3]){
        summaries = fread(file.name)
        summaries = summaries[, c("#CHROM", "POS", "ID", "REF", "ALT", "A1", "A1_FREQ", "OBS_CT", "BETA", "SE")]
        names(summaries)[c(ncol(summaries)-1, ncol(summaries))] = c(paste(i, "BETA", sep="_"), paste(i, "SE", sep="_"))
      } else {
        res = fread(file.name, select = c("BETA", "SE"))
        summaries$newbeta = res$BETA
        summaries$newse = res$SE
        names(summaries)[c(ncol(summaries)-1, ncol(summaries))] = c(paste(i, "BETA", sep="_"), paste(i, "SE", sep="_"))
      }
      print(sprintf("Done for %s", i))
    }
    n_phenos = dim(phenos_of_different_sizes[[j]])[2]-2
    n_individuals = dim(phenos_of_different_sizes[[j]])[1]
    matrix_name = sprintf("%d_individuals_%d_phenos", n_individuals, n_phenos)
    all_summaries[[matrix_name]] = summaries
    matrix_names[[j]] = matrix_name
    
    beta_idxs <- seq(9, dim(summaries)[2]-1, by = 2)
    se_idxs <- seq(10,  dim(summaries)[2], by = 2)
    snp_ids <- summaries$ID
    
    betas <- summaries[, ..beta_idxs]
    se <- summaries[, ..se_idxs]
    
    # Calculate p-values
    chisq_stats <- (betas / se)^2
    p_vals <- pchisq(as.matrix(chisq_stats), df = 1, lower.tail = FALSE)
    
    # Selecting only those which do not have an effect
    idxs_to_keep <- which(colMins(p_vals) > threshold)
    betas_reduced <- betas[,..idxs_to_keep]
    num_traits_kept <- ncol(data.frame(betas_reduced))
    percentage_phenos <- (num_traits_kept/pheno_sizes[j]) * 100
    se_reduced <- se[,..idxs_to_keep]
    
    zscores <- matrix(0, nrow = dim(betas_reduced)[1], dim(betas_reduced)[2])
    zscores <- as.data.frame(zscores)
    
    for (i in 1:ncol(betas_reduced)) {
      zscores[[i]] <- betas_reduced[[i]] / se_reduced[[i]]
    }
    
    # Now we calculate the correlation between the betas and the zscores
    betas_cor <- WGCNA::cor(t(betas_reduced))
    rownames(betas_cor) <- snp_ids
    colnames(betas_cor) <- snp_ids
    zscore_cor <- WGCNA::cor(t(zscores)) 
    rownames(zscore_cor) <- snp_ids
    colnames(zscore_cor) <- snp_ids
    
    r_sim_zscore[[matrix_name]] = zscore_cor
    r_sim_betas[[matrix_name]] = betas_cor
    perc_phenos_remaining_post_thresholding[[matrix_name]] = percentage_phenos
  }
  
  # Assuming matrix_names and ld_matrices are already defined
  output <- list(
    `matrix_names_list` = matrix_names,
    `R_sim_zscore_list` = r_sim_zscore,
    `R_sim_betas_list` = r_sim_betas,
    `sum_stats_list` = all_summaries,
    `perc_phenos_left` = perc_phenos_remaining_post_thresholding
  )
  
  system(sprintf("rm *.phen"))
  if (subsample_individuals == TRUE) {
    system(sprintf("rm %s*", b.file.extension.subsample))
  } else {
    system(sprintf("rm %s*", b.file.extension))
  }
  
  saveRDS(output, sprintf("model_%d_with_R_sim.rds", model.no))
  
}

models <- read_table("~/data/2025_13_05_models.txt")
chromosomes <- models$chr[1:10]
start.points <- models$start[1:10]
end.points <- models$end[1:10]

prepare_thresholded_ld (
     model.no=model_number
    ,null.phenotype.path="/scratch/giulia.desanctis/susie_results_method_2/phenotypes_two_thirds_population/model_%d/20k_pheno_file_model_%d.rds"
    ,effect.phenotype.path = "/scratch/giulia.desanctis/susie_results_method_2/phenotypes_two_thirds_population/model_%d/5k_pheno_file_model_%d_effect_%.2f.rds"
    ,pheno_sizes=c(1000, 2000, 3000, 4000, 5000, 6000, 7000, 8000, 9000, 10000)
    ,bed.internal="/ssu/bsssu/giulia_ld_simulations/2025_05_02_brit_iris_genotype_data/ukbb_all_chrs_grch38_maf0.01_30000_random_unrelated_white_british_alpha_sort_alleles_aligned_nodups_for_sim" 
    ,plink2.binary="/ssu/gassu/software/plink/2.00_20211217/plink2"
    ,plink.binary="/ssu/gassu/software/plink/1.90_20210606/plink"
    ,maf.threshold=0.01
    ,chr=chromosomes[model_number]
    ,start.pos=start.points[model_number]
    ,end.pos=end.points[model_number]				
    ,keep.snps=FALSE
    ,snp.list=NULL
    ,subsample_individuals=TRUE
    ,individuals="~/individual_lists/individual_idxs_internal_10k.txt"
    ,num_individuals=10000
    ,b.file.extension="all_indiv"
    ,b.file.extension.subsample="subsample_indiv"
    ,keep.bfile=FALSE
    ,effect.sizes=c(0.2)
    ,dec.phenos.with.effect=0.1
    ,threshold=1e-6
    )

