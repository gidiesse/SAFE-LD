# This is Nicola's original function modified in the phenotype generation aspect.
# There is a new term called "perc.phenos.with.effect" and this term is used to set
# the percentage of phenotypes that you want to be simulated with an effect. e.g. if we 
# set it 0.1, then 10% of phenotypes will be generated with an effect and the remaining 90% without an 
# effect

# For more information about the general function and parameters refer to:
# fun_nicola_cluster_bed_fix.R

#setwd("/home/giulia.desanctis/simulated_LD_project_giulia_ds/LD_simulation_project_experiments/")

set.seed(42)

library(genio)
library(data.table)
library(susieR)
library(pgenlibr)
library(WGCNA)
library(stringr)
library(pracma)
library(abind)
library(matrixStats)

beta.calc=function(r2,maf) sqrt(r2/(2*maf*(1-maf)))

simulate_dataset=function(
     effect_sizes=c(0) 
    #effect_sizes=c(0,0)
    ,bed.files="/processing_data/shared_datasets/ukbiobank/genotypes/LD_reference/ld_reference_bfiles/ukbb_all_30000_random_unrelated_white_british"
    ,maf.threshold=0.01 
    ,chr=22 
    ,start.pos=32833 
    ,end.pos=33083
    ,keep.bfile=FALSE
    ,plink2.binary="/ssu/gassu/software/plink/2.00_20211217/plink2"
    ,b.file.extension="test_giulia_b"
    ,b.file.extension.subsample="test_giulia_b_subsample"
    ,n.traits=1900
    ,num_individuals=1000
    ,keep.snps=FALSE
    ,snp.list=NULL
    ,perc.phenos.with.effect=0.1
    ,chosen_seed=123
    )
{
  
  # Extract region
  if(keep.bfile==FALSE){
    commando1=paste(plink2.binary,"--bfile",bed.files,"--keep-allele-order","--chr",chr,"--from-kb",start.pos,"--to-kb",end.pos,"--maf",maf.threshold, "--make-bed", "--out", b.file.extension)
    system(commando1)
  }

  genos_full = read_plink(b.file.extension)
  
  # Extract number of desired individuals
  set.seed(chosen_seed)
  individuals <- sample(genos_full$fam$fam,num_individuals,FALSE)
  df <- data.frame(FID = individuals, IID = individuals)
  write.table(df, "individual_idxs.txt", sep = " ", row.names = FALSE, col.names = FALSE, quote = FALSE)
  set.seed(NULL)
  
  # Now we filter for MAF and prepare data for the required number of individuals 
  commando1=paste(plink2.binary,"--bfile",b.file.extension,"--keep","individual_idxs.txt","--maf",maf.threshold,"--make-bed", "--out", b.file.extension.subsample)
  system(commando1)
  
  # Now we read the subsampled file
  genos = read_plink(b.file.extension.subsample)
  
  # Check for NAs
  if (sum(is.na(genos$X))!=0) {
    stop("NAs present!")
  }
  
  freqs=rowMeans(genos$X)/2
  effect.list=list()
  pheno.file=data.frame(matrix(NA,nrow=nrow(genos$fam),ncol=n.traits))
  
  num_phenos_with_effect = ceil(n.traits * perc.phenos.with.effect)
  for(i in 1:n.traits){
    if(keep.snps==FALSE | i==1){
      causal_idx=sample(1:nrow(genos$X),size = length(effect_sizes),replace = F)
    }
    if(!is.null(snp.list)){
      causal_idx=which(genos$bim$id%in%snp.list)
    }
    # Added by me to give names to the columns and rows of genos$X#
    c_names <- genos$fam$fam
    colnames(genos$X) <- c_names
    r_names <- genos$bim$id
    rownames(genos$X) <- r_names
    #End of section added by me#
    causal_snps=row.names(genos$X)[causal_idx]
    betas=beta.calc(r2=effect_sizes,maf=freqs[causal_idx])
    effect.list[[i]]=data.frame(snps=causal_snps,effects=betas,r2=effect_sizes)
    if(length(betas)==1){
      tratto.g=t(genos$X[causal_snps,])*betas
    }else{
      tratto.g=t(genos$X[causal_snps,])%*%betas
    }
    if (between(i, 1, num_phenos_with_effect)) {
      tratto=as.vector(tratto.g)+rnorm(n = ncol(genos$X),sd = sqrt(1-(sum(effect_sizes))))
      pheno.file[,i]=tratto
    } else{
      tratto=rnorm(n = ncol(genos$X),sd = sqrt(1))
      pheno.file[,i]=tratto
    }
  }
  
  # Check that everything has been filled correctly
  if (sum(is.na(pheno.file)) != 0) {
    stop("There are NAs in the phenotypes!")
  }
  
  pheno.file = data.frame( FID=genos$fam$fam,IID=genos$fam$id, pheno.file )
  names(pheno.file)=c("FID","IID",paste0("P",1:n.traits))
  write.table(pheno.file,file = paste0(b.file.extension.subsample,".phen"),row.names=F,quote=F,sep="\t")
  command=paste(plink2.binary,"--bfile",b.file.extension.subsample, "--pheno",paste0(b.file.extension.subsample,".phen")," --glm cols=+a1freq allow-no-covars omit-ref --out",b.file.extension.subsample)
  system(command)
  #### results collector
  summaries=c()
  for ( i in colnames(pheno.file)[-c(1:2)]){
    file.name=paste0(b.file.extension.subsample,".",i,".glm.linear")
    if(i == colnames(pheno.file)[3]){
      summaries=fread(file.name)
      summaries=summaries[,c("#CHROM","POS","ID" ,"REF" ,  "ALT"  ,  "A1","A1_FREQ","OBS_CT","BETA","SE")]
      names(summaries)[c(ncol(summaries)-1,ncol(summaries))]=c(paste(i,"BETA",sep="_"),paste(i,"SE",sep="_"))
    }else{
      res=fread(file.name,select = c("BETA","SE"))
      summaries$newbeta=res$BETA
      summaries$newse=res$SE
      names(summaries)[c(ncol(summaries)-1,ncol(summaries))]=c(paste(i,"BETA",sep="_"),paste(i,"SE",sep="_"))
      
    }
  }
  maf.scale=function(x){
    variance = var(x)
    x=(x-mean(x))/sqrt(variance)
  }
  scaled.matr=apply(t(genos$X),2,  maf.scale)
  cor.matr=(t(scaled.matr)%*%scaled.matr)/nrow(genos$fam)
  diag(cor.matr)=1
  
  output=list(effects.list=effect.list,s_stat=as.data.frame(summaries),R=cor.matr,
              N=nrow(genos$fam), n.traits=n.traits, 
              perc_phenos_with_effect=perc.phenos.with.effect,
              effect_size=effect_sizes, num_individuals=num_individuals)
  return(output)
}

process_dataset=function(
    bed.files, num_it, thresholds, effect.size, n.traits, num.individuals, 
    perc.phenos.with.effect,start.pos=32833,end.pos=33083,chr=22
    
)
  
{
  results <- data.frame(matrix(NA, nrow = length(thresholds), ncol = 2*num_it), 
                        row.names = thresholds)
  for (it in 1:num_it) {
    
    start.time <- Sys.time()
    output = simulate_dataset(bed.files=bed.files,effect_sizes=effect.size,n.traits=n.traits ,num_individuals=num.individuals, perc.phenos.with.effect=perc.phenos.with.effect)
    print(sprintf("Using these files: %s", bed.files))
    cor_mat = output$R
    summaries <- output$s_stat
    
    # Capture the betas and the SE 
    beta_idxs <- seq(9, dim(summaries[,,1])[2]-1, by = 2)
    se_idxs <- seq(10,  dim(summaries[,,1])[2], by = 2)
    snp_ids <- output$s_stat[,3]
    
    betas <- output$s_stat[, beta_idxs]
    se <- output$s_stat[, se_idxs]
    
    write.table(cor_mat, "cor_mat.txt", row.names=TRUE, col.names = TRUE, sep = "\t")
    
    chisq_stats <- (betas / se)^2
    p_vals <- pchisq(as.matrix(chisq_stats), df = 1, lower.tail = FALSE)
    
    for (i in 1:length(thresholds)) {
      idxs_to_keep <- which(colMins(p_vals) > thresholds[i])
      betas_reduced <- betas[,idxs_to_keep]
      num_traits_kept <- ncol(data.frame(betas_reduced))
      percentage_phenos <- (num_traits_kept/n.traits) * 100
      se_reduced <- se[,idxs_to_keep]
      
      if (dim(data.frame(betas_reduced))[2] > 1) {
        # correlate betas
        betas_cor <- WGCNA::cor(t(betas_reduced))
        row.names(betas_cor) <- snp_ids
        colnames(betas_cor) <- snp_ids
        
        # Save this 
        write.table(betas_cor, "cor_betas.txt", row.names=TRUE, col.names = TRUE, sep = "\t")
        write.table(betas_reduced, "betas_reduced.txt", row.names=TRUE, col.names = TRUE, sep = "\t")
        write.table(se_reduced, "se_reduced.txt", row.names=TRUE, col.names = TRUE, sep = "\t")
        
        # Turn matrices into a lower triangular matrix 
        lu_cor_mat <- cor_mat[lower.tri((cor_mat))]
        lu_betas_cor <- betas_cor[lower.tri((betas_cor))]
        
        percentage_phenos_effect = perc.phenos.with.effect * 100
        
        # Do plot of betas against cor 
        jpeg(file="snps_betas_cor.jpg")
        plot(betas_cor,cor_mat)
        abline(a = 0, b = 1, col = "red", lwd = 2)
        title(main = sprintf("Cor(snps) vs cor(betas) - %d phenos %d individuals it %d \n - %f effect on %.2f percent phenos %.10f p-val threshold", 
                             as.integer(n.traits), as.integer(num.individuals), as.integer(it), effect.size,  percentage_phenos_effect, thresholds[i]))
        dev.off()
        
        # Do boxplot
        jpeg(file="box_plot.jpg")
        diff <- betas_cor - cor_mat
        boxplot((as.vector(diff))~findInterval(as.vector(cor_mat),vec=seq(from=-1,to=1,by=0.1)))
        title(main = sprintf("Boxplot - %d phenos %d individuals it %d \n %f effect on %.2f percent phenos %.10f p-val threshold", 
                             as.integer(n.traits),  as.integer(num.individuals), as.integer(it), effect.size,  percentage_phenos_effect, thresholds[i]))
        dev.off()
        
        # save correlation
        cor_mat <- as.vector(cor_mat)
        betas_cor <- as.vector(betas_cor)
        correlation <- cor(cor_mat, betas_cor)
        
        colnames(results)[2*it-1] <- sprintf("Corr. it_%d", it)
        colnames(results)[2*it] <- sprintf("perc phenos left it_%d", it)
        results[i,2*it-1] <- correlation
        results[i,2*it] <- percentage_phenos
      } else {
        #betas_cor <- matrix(ncol = 0, nrow = 0)
        colnames(results)[2*it-1] <- sprintf("Corr. it_%d", it)
        colnames(results)[2*it] <- sprintf("perc phenos left it_%d", it)
        results[i,2*it-1] <- "N/A"
        results[i,2*it] <- percentage_phenos
      }
      system(sprintf("mkdir threshold_%.10f", thresholds[i]))
      system(sprintf("mv snps_betas_cor.jpg threshold_%.10f", thresholds[i]))
      system(sprintf("mv box_plot.jpg threshold_%.10f", thresholds[i]))
      system(sprintf("mv betas_reduced.txt threshold_%.10f", thresholds[i]))
      system(sprintf("mv cor_betas.txt threshold_%.10f", thresholds[i]))
      system(sprintf("mv se_reduced.txt threshold_%.10f", thresholds[i]))
    }
  
    
    if (it == 1) {
      system(sprintf("mkdir %d_iterations_%d_phenos_%f_percent_with_effect_size_of_%f_UKB", 
                     as.integer(num_it), as.integer(n.traits), perc.phenos.with.effect, effect.size))
    }
    
    system(sprintf("mkdir it_%d", it))
    system(sprintf("mv test_giulia_b.log it_%d", it))
    system(sprintf("mv individual_idxs.txt it_%d", it))
    system("rm test_giulia_b_subsample.*")
    system("rm test_giulia_b.*")
    system(sprintf("mv threshold_* it_%d", it))
    
    cor_idxs <- seq(1, 2*num_it, 2)
    pheno_idxs <- seq(2, 2*num_it, 2)
    results$average_cor <- rep(NA, length(thresholds))
    results$average_perc_phenos_left <- rep(NA, length(thresholds))
    results$variance_cor <- rep(NA, length(thresholds))
    for (i in 1:length(thresholds)) {
      cors <- results[i, cor_idxs]
      perc_phenos <- results[i, pheno_idxs]
      results[i,"average_cor"] <- mean(as.numeric(cors), na.rm = TRUE)
      results[i,"average_perc_phenos_left"] <- mean(as.numeric(perc_phenos), na.rm = TRUE)
      if (is.na(var(as.numeric(cors)))) {
        results[i,"variance_cor"] <- "N/A, just one cor or no cor are non-zero"
      } else {
        results[i,"variance_cor"] <- var(as.numeric(cors), na.rm = TRUE)
      }
    }
    write.csv(results, "results.csv", row.names = TRUE)
    
    if (it == num_it ) {
      for (g in 1:num_it) {
        system(sprintf("mv it_%d %d_iterations_%d_phenos_%f_percent_with_effect_size_of_%f_UKB", 
                       as.integer(g), as.integer(num_it), as.integer(n.traits), 
                       perc.phenos.with.effect, effect.size))
      }
      system(sprintf("mv cor_mat.txt %d_iterations_%d_phenos_%f_percent_with_effect_size_of_%f_UKB", 
                     as.integer(num_it), as.integer(n.traits), 
                     perc.phenos.with.effect, effect.size))
      system(sprintf("mv results.csv %d_iterations_%d_phenos_%f_percent_with_effect_size_of_%f_UKB", 
                     as.integer(num_it), as.integer(n.traits), 
                     perc.phenos.with.effect, effect.size))
    }
    
    end.time <- Sys.time()
    total.time <- end.time - start.time
    print(sprintf("Time taken: %f", total.time))
  }
}

process_dataset(bed.files="/ssu/bsssu/ghrc38_reference/ukbb_all_chrs_grch38_maf0.01_30000_random_unrelated_white_british_alpha_sort_alleles", 
                num_it=5, thresholds=c(1e-8,1e-7,1e-6,1e-5,1e-4), effect.size=c(0.01), n.traits=2000, num.individuals=1000, 
                  perc.phenos.with.effect=0.01,start.pos=23175,end.pos=24175,chr=6)
