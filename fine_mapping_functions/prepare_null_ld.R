library(genio)
library(data.table)
library(susieR)
library(pgenlibr)
library(WGCNA)
library(stringr)
library(dplyr)

beta.calc=function(r2,maf) sqrt(r2/(2*maf*(1-maf)))
maf.scale=function(x){
  variance = var(x)
  x=(x-mean(x))/sqrt(variance)
}

prepare_null_ld <- function(
     model.no=2
    ,phenotype.path="/project/statgen/models_for_final_simulations/model_%d/20k_pheno_file_model_%d.rds"
    ,pheno_sizes=c(5000, 3000, 1000, 50, 10)
    ,bed.internal="/ssu/bsssu/giulia_ld_simulations/2025_05_02_brit_iris_genotype_data/ukbb_all_chrs_grch38_maf0.01_30000_random_unrelated_white_british_alpha_sort_alleles_aligned_nodups_for_sim" 
    ,plink2.binary="/ssu/gassu/software/plink/2.00_20211217/plink2"
    ,plink.binary="/ssu/gassu/software/plink/1.90_20210606/plink"
    ,maf.threshold=0.01
    ,chr=19
    ,start.pos=44745	
    ,end.pos=45745			
    ,keep.snps=FALSE
    ,snp.list=NULL
    ,subsample_individuals=TRUE
    ,individuals="/home/giulia.desanctis/fine_mapping/5k_individuals.txt"
    ,num_individuals=5000
    ,b.file.extension="all_indiv"
    ,b.file.extension.subsample="subsample_indiv"
    ,keep.bfile=FALSE
    ,effect.sizes=c(0)
    ) 
{
  
  # 1) Using internal bed files for 30k individuals
  
  all_phenos <- readRDS(sprintf(phenotype.path, model.no, model.no))
  
  # Now we need to subset the phenotypes
  phenos_of_different_sizes <- list()
  
  for (i in 1:length(pheno_sizes)) {
    idxs <- sample(3:20002, pheno_sizes[i], replace = FALSE)
    phenos <- all_phenos$phenos[,c(1,2,idxs)]
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
    
    #scaled.matr=apply(t(genos_full$X),2,  maf.scale)
    #cor.matr=(t(scaled.matr)%*%scaled.matr)/nrow(genos_full$fam)
    #diag(cor.matr)=1
  
    #ld_matrices[[matrix_name]] = cor.matr
    
    beta_idxs <- seq(9, dim(summaries)[2]-1, by = 2)
    se_idxs <- seq(10,  dim(summaries)[2], by = 2)
    snp_ids <- summaries$ID
    
    betas <- summaries[, ..beta_idxs]
    se <- summaries[, ..se_idxs]
    zscores <- matrix(0, nrow = dim(betas)[1], dim(betas)[2])
    zscores <- as.data.frame(zscores)
    
    for (i in 1:pheno_sizes[j]) {
      zscores[[i]] <- betas[[i]] / se[[i]]
    }
    
    # Now we calculate the correlation between the betas and the zscores
    betas_cor <- WGCNA::cor(t(betas))
    rownames(betas_cor) <- snp_ids
    colnames(betas_cor) <- snp_ids
    zscore_cor <- WGCNA::cor(t(zscores)) 
    rownames(zscore_cor) <- snp_ids
    colnames(zscore_cor) <- snp_ids
    
    r_sim_zscore[[matrix_name]] = zscore_cor
    r_sim_betas[[matrix_name]] = betas_cor
  }
  
  # Assuming matrix_names and ld_matrices are already defined
  output <- list(
    `matrix_names_list` = matrix_names,
    `R_sim_zscore_list` = r_sim_zscore,
    `R_sim_betas_list` = r_sim_betas,
    `sum_stats_list` = all_summaries
  )
  
  system(sprintf("rm *.phen"))
  system(sprintf("rm %s*", b.file.extension))
  system(sprintf("rm %s*", b.file.extension.subsample))
  
  saveRDS(output, sprintf("model_%d_with_R_sim_more_phenos_than_indiv.rds", model.no))
  
}

model_info <- fread("~/data/2025_13_05_models.txt")  
chromosomes <- model_info$chr[1:100]
starts <- model_info$start[1:100]
ends <- model_info$end[1:100]
   
for (i in 1:10) {
  
  setwd(sprintf("/scratch/giulia.desanctis/susie_results/5000_individuals/model_%d", i))
  
  start.time <- Sys.time()
  prepare_null_ld (
     model.no=i
    ,phenotype.path="/project/statgen/models_for_final_simulations/model_%d/20k_pheno_file_model_%d.rds"
    ,pheno_sizes=c(5000, 6000, 7000, 8000, 9000, 10000)
    ,bed.internal="/scratch/giulia.desanctis/british_data_split_for_int_ext/internal_two_thirds"
    ,plink2.binary="/ssu/gassu/software/plink/2.00_20211217/plink2"
    ,maf.threshold=0.01
    ,chr=chromosomes[i]
    ,start.pos=starts[i]
    ,end.pos=ends[i]
    ,keep.snps=FALSE
    ,snp.list=NULL
    ,subsample_individuals=TRUE
    ,individuals="~/individual_lists/individual_idxs_internal_5k.txt"
    ,num_individuals=5000
    ,b.file.extension="all_indiv"
    ,b.file.extension.subsample="subsample_indiv"
    ,keep.bfile=FALSE
    ,effect.sizes=c(0)
  ) 
  
  end.time <- Sys.time()
  
  print(sprintf("CREATED ALL THE R_SIM MATRICES FOR MODEL %d !!!!", i))
  time.taken <- end.time - start.time
  print(sprintf("time taken: %f", time.taken))

  
}
