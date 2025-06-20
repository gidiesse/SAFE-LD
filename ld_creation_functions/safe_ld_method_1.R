# The OG function and the function that is at the base of the entire project.
# This function will have a very long explanation and then every other function 
# that begins with fun_Nicola"...".R all derive from this function and will be 
# modified for a specific purpose which will be stated at the top of the function.

# This function was written by Nicola and then modified by me to fit the project 
# and it has the following parameters: 
# 1) "effect_sizes" - by default this is at c(0), but to add an effect on the phenotypes
# generated you just put the percentage of variance explained by that phenotype here 
# e.g. c(0.01,0.05) are two snps that explain respectively 1% and 5% of variance 
# 2) "bed.files" - filepath of the bed files 
# 3) "maf.threshold" - what you want the maf to be thresholded at 
# 4) "chr" - when using the --from and --to flags in plink, the chromosome number 
# has to be specified 
# "start.pos" - this is for the --from flag in plink and denotes the starting position of 
# the extracted region 
# "end.pos" - this is for the --to flag in plink and denotes the end position of 
# the extracted region 
# "plink2.binary" - filepath where plink2 is saved 
# "b.file.extension" - name for how the bed files will be saved once they desired
# region is extracted and they're filtered for MAF and HWE equilibrium
# "n.traits" - number of phenotypes to be simulated 
# "keep.snps" - this is for the functioning of the phenotype generation 
# "snp.list" - this can be used if you want specific SNPs to be simulated as causal,
# otherwise it can simply be set to NULL 
# "snp.range" - this is an alternative to extracting a specific region, instead of
# giving start and end position you can use a snp.range to extract a section of the genome
# "seed" - if you want to set a seed for reproducibility of results 

# There are two other much smaller functions defined in this script:
# 1) beta.calc - written by Nicola to calculate the beta coefficients
# 2) mean_of_col_with_NA - written by me when we were dealing with pgen data 
# and there were NAs, this could probably be archived now 

setwd("/project/statgen/results_nulls_different_chromosomes")

library(genio)
library(data.table)
library(susieR)
library(pgenlibr)
library(WGCNA)
library(stringr)
library(filesstrings)
library(fs)

beta.calc=function(r2,maf) sqrt(r2/(2*maf*(1-maf)))

simulate_dataset=function(
     effect_sizes=c(0)
    ,bed.files="/processing_data/shared_datasets/ukbiobank/genotypes/LD_reference/ld_reference_bfiles/ukbb_all_30000_random_unrelated_white_british"
    ,maf.threshold=0.01
    ,chr=1
    ,start.pos=32833   
    ,end.pos=33833 
    ,plink2.binary="/ssu/gassu/software/plink/2.00_20211217/plink2"
    ,b.file.extension="test_giulia_b"
    ,n.traits=5000
    ,keep.snps=FALSE
    ,snp.list=NULL
    )
{
  
  commando=paste(plink2.binary,"--bfile",bed.files,"--keep-allele-order","--chr",chr,"--from-kb",start.pos,"--to-kb",end.pos,"--maf",maf.threshold, "--make-bed", "--out", b.file.extension)
  system(commando)

  #commando1=paste(plink2.binary,"--bfile",bed.files,"--keep-allele-order","--snps",snp.range,"--maf",maf.threshold, "--make-bed", "--out", b.file.extension)
  #system(commando1)

  genos = read_plink(b.file.extension)
  
  # Check for NAs
  if (sum(which(is.na(genos$X))) != 0) {
    stop("NAs are present in the genetic data!")
  }
  
  freqs=rowMeans(genos$X)/2
  effect.list=list()
  pheno.file=data.frame(matrix(NA,nrow=nrow(genos$fam),ncol=n.traits))
  for(i in 1:n.traits){
    if(keep.snps==FALSE | i==1){
      causal_idx=sample(1:nrow(genos$X),size = length(effect_sizes),replace = F)
    }
    if(!is.null(snp.list)){
      causal_idx=which(genos$bim$id%in%snp.list)
    }
    
    causal_snps=row.names(genos$X)[causal_idx]
    betas=beta.calc(r2=effect_sizes,maf=freqs[causal_idx])
    effect.list[[i]]=data.frame(snps=causal_snps,effects=betas,r2=effect_sizes)
    if(length(betas)==1){
      tratto.g=t(genos$X[causal_snps,])*betas
    }else{
      tratto.g=t(genos$X[causal_snps,])%*%betas
    }
    tratto=as.vector(tratto.g)+rnorm(n = ncol(genos$X),sd = sqrt(1-(sum(effect_sizes))))
    pheno.file[,i]=tratto
  }
  pheno.file = data.frame( FID=genos$fam$fam,IID=genos$fam$id, pheno.file )
  names(pheno.file)=c("FID","IID",paste0("P",1:n.traits))
  write.table(pheno.file,file = paste0(b.file.extension,".phen"),row.names=F,quote=F,sep="\t")
  command=paste(plink2.binary,"--bfile",b.file.extension, "--pheno",paste0(b.file.extension,".phen")," --glm cols=+a1freq allow-no-covars omit-ref --out",b.file.extension)
  system(command)
  #### results collector
  start <- Sys.time()
  # Initialize summary table
  summaries = NULL
  colnames_to_iterate = colnames(pheno.file)[-c(1:2)]
  
  # Process first column separately to initialize 'summaries'
  first_col = colnames_to_iterate[1]
  file.name = paste0(b.file.extension, ".", first_col, ".glm.linear")
  
  summaries = fread(file.name, select = c("#CHROM", "POS", "ID", "REF", "ALT", "A1", "A1_FREQ", "OBS_CT", "BETA", "SE"))
  setnames(summaries, c("BETA", "SE"), c(paste0(first_col, "_BETA"), paste0(first_col, "_SE")))
  
  # Process remaining columns efficiently
  for (i in colnames_to_iterate[-1]) {
    file.name = paste0(b.file.extension, ".", i, ".glm.linear")
    res = fread(file.name, select = c("BETA", "SE"))
    summaries[, paste0(i, "_BETA") := res$BETA]
    summaries[, paste0(i, "_SE") := res$SE]
  }
  
  # Scaling function using vectorized operations
  maf.scale <- function(x) {
    (x - mean(x, na.rm = TRUE)) / sqrt(var(x, na.rm = TRUE))
  }
  
  # Apply scaling efficiently
  scaled.matr = scale(t(genos$X))  # More efficient than apply()
  cor.matr = crossprod(scaled.matr) / nrow(genos$fam)
  diag(cor.matr) = 1
  
  # Output
  output = list(effects.list = effect.list, s_stat = summaries, R = cor.matr, N = nrow(genos$fam))
  
  end= Sys.time()
  time_old_2 = end - start
  print("time taken for results collection:")
  print(time_old_2)
  return(output)
}


# To run multiple iterations 
num_it = 100

correlation_values <- vector(mode = "numeric", length = num_it)

for (k in 1:num_it) {
  
  start.time <- Sys.time()

  output = simulate_dataset()
  
  end.time <- Sys.time()
  execution_time <- end.time - start.time
  print(execution_time)
  
  cor_mat = output$R
  n.traits = length(output$effects.list)
  
  # Now we want to capture the betas and the SE so we can make plots 
  
  beta_idxs <- seq(9, dim(output$s_stat)[2]-1, by = 2)
  se_idxs <- seq(10,  dim(output$s_stat)[2], by = 2)
  snp_ids <- output$s_stat[,3]
  
  write.table(cor_mat, "cor_mat.txt", row.names=TRUE, col.names = TRUE, sep = "\t")
  
  betas <- as.matrix(output$s_stat[, beta_idxs, with = FALSE])
  se <- as.matrix(output$s_stat[, se_idxs, with = FALSE])
  zscores <- betas / se
  
  
  # Now we calculate the correlation between the betas and the zscores
  betas_cor <- WGCNA::cor(t(betas))
  row.names(betas_cor) <- snp_ids[[1]]
  colnames(betas_cor) <- snp_ids[[1]]
  zscore_cor <- WGCNA::cor(t(zscores)) 
  row.names(zscore_cor) <- snp_ids[[1]]
  colnames(zscore_cor) <- snp_ids[[1]]
  
  # Save this 
  write.table(betas_cor, "cor_betas.txt", row.names=TRUE, col.names = TRUE, sep = "\t")
  write.table(zscore_cor, "cor_zscore.txt", row.names=TRUE, col.names = TRUE, sep = "\t")
  
  # Calculate correlation
  zscore_cor_mat <- as.matrix(zscore_cor)
  zscore_cor_vec <- as.vector(zscore_cor_mat)
  cor_mat_mat <- as.vector(cor_mat)
  cor_mat_vec <- as.vector(cor_mat_mat)
  
  correlation_values[k] <- cor(cor_mat_vec, zscore_cor_vec)
  
  if (k == 1) {
  # Turn matrices into a lower triangular matrix 
  lu_cor_mat <- cor_mat[lower.tri((cor_mat))]
  
  lu_betas_cor <- betas_cor[lower.tri((betas_cor))]
  lu_zscore_cor <- zscore_cor[lower.tri((zscore_cor))]
  
  # Do plot of betas against cor 
  jpeg(file="snps_betas_cor.jpg")
  plot(betas_cor,cor_mat)
  title(sprintf("Cor(snps) vs cor(betas) UKB - %d phenos - it %d", n.traits, k))
  dev.off()
  
  # Do plot of zscore against cor 
  jpeg(file="snps_zscore_cor.jpg")
  plot(zscore_cor,cor_mat)
  title(sprintf("Cor(snps) vs cor(zscore) UKB - %d phenos - it %d", n.traits, k))
  dev.off()
  
  # Do boxplot
  jpeg(file="box_plot.jpg")
  diff = betas_cor - cor_mat
  boxplot((as.vector(diff))~findInterval(as.vector(output$R),vec=seq(from=-1,to=1,by=0.1)))
  title(sprintf("Boxplot - %d phenos - it %d - nulls", n.traits, k))
  dev.off()
  }
  
  dir.create(sprintf("it_%d", k))
  
  if (k == 1) {
    dir.create(sprintf("%d_iterations_%d_phenos_ukb_3k_phenos_basic_script", num_it, n.traits))
    move_files("cor_mat.txt", sprintf("it_%d", k))
    move_files("snps_betas_cor.jpg", sprintf("it_%d", k))
    move_files("snps_zscore_cor.jpg", sprintf("it_%d", k))
    move_files("box_plot.jpg", sprintf("it_%d", k))
  }
 
  move_files("cor_betas.txt", sprintf("it_%d", k))
  move_files("cor_zscore.txt", sprintf("it_%d", k))
  move_files("test_giulia_b.log", sprintf("it_%d", k))
  remove_dir("test_giulia_b.*")
  remove_dir("cor_mat.txt")
  
  if (k == num_it ) {
    for (g in 1:num_it) {
      fs::file_move(sprintf("it_%d", g), sprintf("%d_iterations_%d_phenos_ukb_3k_phenos_basic_script", num_it, n.traits))
      
    }
    write(correlation_values, file = "correlation_real_LD_sim_LD_scores.txt")
    move_files("correlation_real_LD_sim_LD_scores.txt", sprintf("%d_iterations_%d_phenos_ukb_3k_phenos_basic_script", num_it, n.traits))
    
  }
  
  end.time <- Sys.time()
  execution_time <- end.time - start.time
  print("time for one iteration: ")
  print(execution_time)
}

