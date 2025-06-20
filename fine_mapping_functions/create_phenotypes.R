args <- commandArgs(trailingOnly = TRUE)
model_number <- as.integer(args[1])

library(genio)
library(data.table)
library(susieR)
library(pgenlibr)
library(WGCNA)
library(stringr)
library(dplyr)
library(readr)
library(pracma)

beta.calc=function(r2,maf) sqrt(r2/(2*maf*(1-maf)))


create_phenotypes=function(
     effect_sizes=c(0)
    ,bed.files="/ssu/bsssu/giulia_ld_simulations/2025_05_02_brit_iris_genotype_data/ukbb_all_chrs_grch38_maf0.01_30000_random_unrelated_white_british_alpha_sort_alleles_aligned_nodups_for_sim"
    ,maf.threshold=0.01
    ,chr=17
    ,start.pos=15221816	
    ,end.pos=16221816	 
    ,keep.bfile=FALSE
    ,plink2.binary="/ssu/gassu/software/plink/2.00_20211217/plink2"
    ,b.file.extension="all_indiv"
    ,keep.snps=FALSE
    ,snp.list = NULL
    ,n.traits=20000
    ,model.no=1
)
{
  
  setwd(sprintf("/scratch/giulia.desanctis/susie_results_method_2/phenotypes_two_thirds_population/model_%d", model.no))
  
  #First we extract the 250kb region 
  if(keep.bfile==FALSE){
    commando1=paste(plink2.binary,"--bfile",bed.files,"--keep-allele-order","--chr",chr,"--from-kb",start.pos,"--to-kb",end.pos,"--maf",maf.threshold,"--make-bed", "--out", b.file.extension)
    system(commando1)
  }
  
  genos = read_plink(b.file.extension)
  
  if (sum(is.na(genos$X))!=0) {
    stop("NAs present!")
  }
  
  freqs=rowMeans(genos$X)/2
  effect.list=list()
  pheno.file=data.frame(matrix(NA,nrow=as.numeric(nrow(genos$fam)),ncol=as.numeric(n.traits)))
  
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
    pheno.file[,i]=as.numeric(tratto)
  }
  pheno.file = data.frame( FID=genos$fam$fam,IID=genos$fam$id, pheno.file )
  names(pheno.file)=c("FID","IID",paste0("P",1:n.traits))
  
  output=list(effects.list=effect.list,phenos=pheno.file)
 
  if(effect_sizes == 0) {
    saveRDS(output, file = sprintf("5k_pheno_file_model_%d.rds", model.no))
  } else {
    saveRDS(output, file = sprintf("5k_pheno_file_model_%d_effect_%.2f.rds", model.no, effect_sizes))
  }
  
  system(sprintf("rm %s*", b.file.extension))
  
  print("20k phenotypes saved")
}

models <- read_table("~/data/2025_13_05_models.txt")
chromosomes <- models$chr[1:10]
start.points <- models$start[1:10]
end.points <- models$end[1:10]

create_phenotypes(
   model.no = model_number
  ,effect_sizes=c(0.2)
  ,bed.files="/scratch/giulia.desanctis/british_data_split_for_int_ext/internal_two_thirds"
  ,maf.threshold=0.01
  ,chr=chromosomes[model_number]
  ,start.pos=start.points[model_number]
  ,end.pos=end.points[model_number]			
  ,keep.bfile=FALSE
  ,plink2.binary="/ssu/gassu/software/plink/2.00_20211217/plink2"
  ,b.file.extension="all_indiv"
  ,keep.snps=FALSE
  ,snp.list = NULL
  ,n.traits=5000
)
  
 

