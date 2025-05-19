library(data.table)
library(susieR)
library(pgenlibr)
library(WGCNA)
library(stringr)
library(dplyr)
library(readr)
library(genio)

beta.calc=function(r2,maf) sqrt(r2/(2*maf*(1-maf)))
maf.scale=function(x){
  variance = var(x)
  x=(x-mean(x))/sqrt(variance)
}

#effect_size_list <- list(c(0.005,0.05,0.2), c(0.005, 0.05, 0), c(0.005,0.2,0), c(0.05, 0.2, 0), c(0.005, 0, 0), c(0.05, 0, 0), c(0.2, 0, 0))

external_and_internal_ld <- function(
     model.no=2
    ,phenotype.path="/project/statgen/models_for_final_simulations/model_%d/20k_pheno_file_model_%d.rds"
    ,n.traits = 100
    ,bed.internal="/ssu/bsssu/giulia_ld_simulations/2025_05_02_brit_iris_genotype_data/ukbb_all_chrs_grch38_maf0.01_30000_random_unrelated_white_british_alpha_sort_alleles_aligned_nodups_for_sim"
    ,bed.external="/ssu/bsssu/giulia_ld_simulations/2025_05_02_brit_iris_genotype_data/ukbb_all_chrs_grch38_unrelated_white_irish_alpha_sorted_alleles_aligned_nodups_for_sim"
    ,plink2.binary="/ssu/gassu/software/plink/2.00_20211217/plink2"
    ,plink.binary="/ssu/gassu/software/plink/1.90_20210606/plink"
    ,effect.size=effect_size
    ,maf.threshold=0.01
    ,chr=19
    ,start.pos=44745	
    ,end.pos=45745
    ,keep.snps=FALSE
    ,snp.list=NULL
    ,subsample_individuals=TRUE
    ,keep.bfile=FALSE
    ,b.file.extension.int="all_indiv_int"
    ,b.file.extension.ext="all_indiv_ext"
) 
{
  setwd(sprintf("/project/statgen/models_for_final_simulations/model_%d", model.no))

     # Internal LD 
    if(keep.bfile==FALSE){
      commando=paste(plink2.binary,"--bfile",bed.internal,"--keep-allele-order","--chr",chr,"--from-kb",start.pos,"--to-kb",end.pos,"--maf",maf.threshold,"--out ", b.file.extension.int," --make-bed")
      system(commando)
    }
    
    genos = read_plink(b.file.extension.int)
    freqs = rowMeans(genos$X)/2
    
    effect.list=list()
    pheno.file=data.frame(matrix(NA,nrow=nrow(genos$fam),ncol=n.traits))
    
    for(i in 1:n.traits){
      
      if(keep.snps==FALSE | i==1){
        causal_idx=sample(1:nrow(genos$X),size = length(effect.size),replace = F)
      }
      if(!any(is.null(snp.list))){
        causal_idx=which(genos$bim$id%in%snp.list)
        
      }
      
      causal_snps=row.names(genos$X)[causal_idx]
      betas=beta.calc(r2=abs(effect.size),maf=freqs[causal_idx]) * sign(effect.size)
      effect.list[[i]]=data.frame(snps=causal_snps,effects=betas,r2=effect.size)
      if(length(betas)==1){
        tratto.g=t(genos$X[causal_snps,])*betas
        
      }else{
        tratto.g=t(genos$X[causal_snps,])%*%betas
      }
      tratto=as.vector(tratto.g)+rnorm(n = ncol(genos$X),sd = sqrt(1-(sum(effect.size))))
      pheno.file[,i]=tratto
    }
    
    pheno.file = data.frame( FID=genos$fam$fam,IID=genos$fam$id, pheno.file )
    names(pheno.file)=c("FID","IID",paste0("P",1:n.traits))
    write.table(pheno.file,file = paste0(b.file.extension.int,".phen"),row.names=F,quote=F,sep="\t")
    command=paste(plink2.binary,"--bfile",b.file.extension.int,"--pheno",paste0(b.file.extension.int,".phen")," --glm cols=+a1freq allow-no-covars omit-ref --out",b.file.extension.int)
    system(command)
    print("association complete")
    
    summaries=c()
    for ( i in colnames(pheno.file)[-c(1:2)]){
      
      file.name=paste0(b.file.extension.int,".",i,".glm.linear")
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
      print(sprintf("summary for %s done", i))
      
    }
    
    summaries <- summaries %>% filter(!duplicated(ID))
    
    maf.scale=function(x){
      
      variance=var(x)
      x=(x-mean(x))/sqrt(variance)
      
    }
    
    scaled.matr.int=apply(t(genos$X),2,  maf.scale)
    cor.matr.int=(t(scaled.matr.int)%*%scaled.matr.int)/nrow(genos$fam)
    diag(cor.matr.int)=1
    
    # External LD
    
    if(keep.bfile==FALSE){
      commando=paste(plink.binary,"--bfile",bed.external,"--chr",chr,"--from-kb",start.pos,"--to-kb",end.pos,"--maf",maf.threshold,"--keep-allele-order", "--out ", b.file.extension.ext," --make-bed")
      system(commando)
    }
    
    genos_ext = read_plink(b.file.extension.ext)
    
    scaled.matr.ext = apply(t(genos_ext$X),2,  maf.scale)
    cor.matr.ext =(t(scaled.matr.ext)%*%scaled.matr.ext)/nrow(genos_ext$fam)
    diag(cor.matr.ext)=1
    
    snps_both <- intersect(row.names(cor.matr.int), row.names(cor.matr.ext))
    cor.matr.int <- cor.matr.int[snps_both, snps_both]
    cor.matr.ext <- cor.matr.ext[snps_both, snps_both]
    
    summaries <- summaries[summaries$ID %in% snps_both, ]
  
  system("rm all_indiv*")
  
  output <- list(
    `R_internal` = cor.matr.int,
    `R_external` = cor.matr.ext,
    `effect_list` = effect.list,
    `sum_stats` = summaries,
    `N_int` = dim(genos$X)[2],
    `N_ext` = dim(genos_ext$X)[2]
  )
    
  saveRDS(output, sprintf("model_%d_with_R_int_R_ext_and_effects.rds", model.no))
  
}

models <- read_table("~/data/2025_13_05_models.txt")
chromosomes <- models$chr
start.points <- models$start
end.points <- models$end
effect.sizes <- models$effect_sizes


for(i in 601:800) {
  
  setwd(sprintf("/project/statgen/models_for_final_simulations/model_%d", i))
  
  external_and_internal_ld(
     model.no=i
    ,phenotype.path="/project/statgen/models_for_final_simulations/model_%d/20k_pheno_file_model_%d.rds"
    ,n.traits = 100
    ,bed.internal="/ssu/bsssu/giulia_ld_simulations/2025_05_02_brit_iris_genotype_data/ukbb_all_chrs_grch38_maf0.01_30000_random_unrelated_white_british_alpha_sort_alleles_aligned_nodups_for_sim"
    ,bed.external="/ssu/bsssu/giulia_ld_simulations/2025_05_02_brit_iris_genotype_data/ukbb_all_chrs_grch38_unrelated_white_irish_alpha_sorted_alleles_aligned_nodups_for_sim"
    ,plink2.binary="/ssu/gassu/software/plink/2.00_20211217/plink2"
    ,plink.binary="/ssu/gassu/software/plink/1.90_20210606/plink"
    ,effect.size= as.numeric(strsplit(effect.sizes, ",")[[i]])
    ,maf.threshold=0.01
    ,chr=chromosomes[i]
    ,start.pos=start.points[i]	
    ,end.pos=end.points[i]
    ,keep.snps=FALSE
    ,snp.list=NULL
    ,subsample_individuals=TRUE
    ,keep.bfile=FALSE
    ,b.file.extension.int="all_indiv_int"
    ,b.file.extension.ext="all_indiv_ext"
  ) 
  
  print(sprintf((" model %d created"),i))
}

