library(genio)
library(data.table)
library(susieR)
library(pgenlibr)
library(WGCNA)
library(stringr)
library(dplyr)
library(readr)

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
    ,b.file.extension.subsample="subsample_indiv"
    ,keep.snps=FALSE
    ,snp.list = NULL
    ,n.traits=20000
    ,model.no=1
)
{
  #First we extract the 250kb region 
  if(keep.bfile==FALSE){
    commando1=paste(plink2.binary,"--bfile",bed.files,"--keep-allele-order","--chr",chr,"--from-kb",start.pos,"--to-kb",end.pos,"--maf",maf.threshold,"--make-bed", "--out", b.file.extension)
    system(commando1)
  }
  
  # Then we extract the number of individuals that we want
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
 
  saveRDS(output, file = sprintf("20k_pheno_file_model_%d.rds", model.no))
  
  system(sprintf("rm %s*", b.file.extension))
  
  print("20k phenotypes saved")
}

models <- read_table("~/data/2025_12_05_models.txt")
chromosomes <- models$chr[1:100]
start.points <- models$start[1:100]
end.points <- models$end[1:100]


for(i in 1:100) {
  
  setwd(sprintf("/project/statgen/models_for_final_simulations/model_%d", i))
  
  create_phenotypes(
     effect_sizes=c(0)
    ,bed.files="/ssu/bsssu/giulia_ld_simulations/2025_05_02_brit_iris_genotype_data/ukbb_all_chrs_grch38_maf0.01_30000_random_unrelated_white_british_alpha_sort_alleles_aligned_nodups_for_sim"
    ,maf.threshold=0.01
    ,chr=chromosomes[i]
    ,start.pos=start.points[i]
    ,end.pos=end.points[i]			
    ,keep.bfile=FALSE
    ,plink2.binary="/ssu/gassu/software/plink/2.00_20211217/plink2"
    ,b.file.extension="all_indiv"
    ,b.file.extension.subsample="subsample_indiv"
    ,keep.snps=FALSE
    ,snp.list = NULL
    ,n.traits=20000
    ,model.no=i
)
  
  print(sprintf(("20 k phenotypes done for model %d"),i))
}

