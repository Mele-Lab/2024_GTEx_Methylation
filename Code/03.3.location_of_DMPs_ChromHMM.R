#!/usr/bin/env Rscript
# @Author: Winona Oliveros Diez
# @E-mail: winn95@gmail.com
# @Description: Get fisher enrichments of DMPs per chromatin state and tissue
# @software version: R=4.2.2


###### Do enrichment of DMPs for each chromatin state per tissue and trait #####

library(dplyr)
library(tidyr)
library(tidyverse)


first_dir <- "/gpfs/projects/bsc83/"
annotation <- read.csv("/gpfs/projects/bsc83/Projects/GTEx_v8/Methylation/Data/GPL21145_MethylationEPIC_15073387_v-1-0_processed.csv")

project_path <- paste0(first_dir, "Projects/GTEx_v8/Methylation/")

tissues <- c("Lung", "ColonTransverse", "Ovary", "Prostate", "BreastMammaryTissue", "MuscleSkeletal", "KidneyCortex", "Testis", "WholeBlood")
names <- c("Age", "Ancestry", "BMI", "Sex")
traits_to_use <- c('EURv1','SEX2','AGE','BMI')

tissues <- c('Ovary')

results_DML <- lapply(tissues, function(tis) 
  readRDS(paste0("/gpfs/projects/bsc83/Projects/GTEx_v8/Methylation/Tissues/",tis,"/DML_results_5_PEERs_continous.rds")))
names(results_DML) <- tissues

files <- list.files('/gpfs/projects/bsc83/Projects/GTEx_v8/Methylation/Data/EpiMap/', pattern='.bed.gz',full.names=T)

names_chrom <- c("Lung", "ColonTransverse", "Ovary", "Prostate", "Breast", "MuscleSkeletal", "KidneyCortex", "Testis", "PBMC")
names_chrom <- c('Ovary')
chromhmm <- lapply(names_chrom, function(tis)
  read.delim(files[grep(tis, files)], sep='\t', header=F))
names(chromhmm) <- tissues

ann_bed <- annotation[!is.na(annotation$MAPINFO) & !is.na(annotation$CHR),] %>%
  dplyr::select(chrom=CHR, start=MAPINFO, end=MAPINFO, name=IlmnID) %>% distinct()
ann_bed$chrom <- paste0('chr',ann_bed$chrom)
head(ann_bed)
ann_bed$start <- ann_bed$start-1

library(valr)
chromhmm_cpgs <- lapply(tissues, function(tis) {
  chrom_df <- chromhmm[[tis]][,c(1:4)]
  colnames(chrom_df) <- c('chrom','start','end','region')
  bed_intersect(ann_bed, chrom_df, suffix = c("_ann", "_chromhmm"))})
names(chromhmm_cpgs) <- tissues

#### enrichment -- chi-squared first ####
### read sharing ####
sharing <- readRDS('/gpfs/projects/bsc83/Projects/GTEx_v8/Methylation/Data/Sharing_DMP.rds')
shared_cpgs <- readRDS('/gpfs/projects/bsc83/Projects/GTEx_v8/Methylation/Data/chromHMM_shared_9tissues.rds')

#From 18 states to 14

# enrichment of DSE in particular types of AS events
my_fisher <- function(type, tissue, trait){
  #             DS      Not DS
  # type
  # other_types
  print(type)
  print(tissue)
  res <- results_DML[[tissue]][[trait]]
  chrom_tissue <- as.data.frame(chromhmm_cpgs[[tissue]])
  chrom_tissue$region_chromhmm_new <- chrom_tissue$region_chromhmm
  chrom_tissue <- chrom_tissue[chrom_tissue$.overlap ==1,]
  # 
  chrom_tissue$region_chromhmm_new[chrom_tissue$region_chromhmm %in% c("TssFlnkD", "TssFlnk", "TssFlnkU","TssA")] <- "TSS"
  chrom_tissue$region_chromhmm_new[chrom_tissue$region_chromhmm %in% c("EnhA2", "EnhA1","EnhWk","EnhG1", "EnhG2")] <- "Enh"

  chrom_tissue$region_chromhmm_new[chrom_tissue$region_chromhmm %in% c("ReprPCWk","ReprPC")] <- "ReprPC"

  chrom_tissue$region_chromhmm_new[chrom_tissue$region_chromhmm %in% c("TxWk","Tx")] <- "Tx"
  
  type_df <- chrom_tissue[chrom_tissue$region_chromhmm_new == type & chrom_tissue$name_ann %in% rownames(res),]
  other_type <- chrom_tissue[chrom_tissue$region_chromhmm_new != type & chrom_tissue$name_ann %in% rownames(res),]
  type_diff <- nrow(type_df[type_df$name_ann %in% rownames(res[res$adj.P.Val<0.05 & res$logFC<0,]),]) # change for hyper
  type_notdiff <- nrow(type_df) - type_diff
  other_type_diff <- nrow(other_type[other_type$name_ann %in% rownames(res[res$adj.P.Val<0.05 & res$logFC<0,]),]) # change for hyper
  other_type_notdiff <- nrow(other_type) - other_type_diff
  
  ### test significance
  m <- matrix(c(type_diff, type_notdiff, other_type_diff, other_type_notdiff), 2,2, byrow = T)
  print(m)  
  
  m[is.na(m)] <- 0
  #m <- m[c(type,paste0('No ',type)),]
  rownames(m) <- c(type, "Other")
  colnames(m) <- c("Hyper","Not Hyper")
  print(m)
  f <- fisher.test(m)
  print(f)
  return(list("f" = f, "m" = type_diff))
}
#Two-tailed Fisher test
families <- as.vector(unique(chromhmm_cpgs$Lung$region_chromhmm))
families <- c('Enh','EnhBiv','Het','Quies','ReprPC','TSS','TssBiv','Tx','ZNF/Rpts')
fisher_results <- lapply(tissues, function(tis) lapply(names(results_DML[[tis]]), function(trait) lapply(families, function(region) my_fisher(region, tis,trait))))
names(fisher_results) <- tissues

for (name in tissues) {
  names(fisher_results[[name]]) <- names(results_DML[[name]])
}

for (name in tissues) {
  for (trait in names(fisher_results[[name]])) {
    names(fisher_results[[name]][[trait]]) <- families
  }
}

saveRDS(fisher_results, '/gpfs/projects/bsc83/Projects/GTEx_v8/Methylation/Tissues/Ovary.enrichment_chromhmm_hypo_batch_CI.simple.continous.filt.rds')

### enrichment shared positions
my_fisher <- function(type, trait){
  #             DS      Not DS
  # type
  # other_types
  print(type)
  res <- sharing[sharing$trait == trait,]
  chrom_tissue <- shared_cpgs

  type_df <- chrom_tissue[chrom_tissue$region_chromhmm == type & chrom_tissue$name_ann %in% res$CG,]
  other_type <- chrom_tissue[chrom_tissue$region_chromhmm != type & chrom_tissue$name_ann %in% res$CG,]
  type_diff <- nrow(type_df[type_df$name_ann %in% res$CG[res$number>=2 & res$dir==1],]) # change for hypo
  type_notdiff <- nrow(type_df) - type_diff
  other_type_diff <- nrow(other_type[other_type$name_ann %in% res$CG[res$number>=2 & res$dir==1],]) # change for hypo
  other_type_notdiff <- nrow(other_type) - other_type_diff

  ### test significance
  m <- matrix(c(type_diff, type_notdiff, other_type_diff, other_type_notdiff), 2,2, byrow = T)
  print(m)

  m[is.na(m)] <- 0
  #m <- m[c(type,paste0('No ',type)),]
  rownames(m) <- c(type, "Other")
  colnames(m) <- c("Hyper","Not Hyper")
  print(m)
  f <- fisher.test(m)
  print(f)
  return(list("f" = f, "m" = type_diff))
}
# Two-tailed Fisher test
#families <- as.vector(unique(shared_cpgs$region_chromhmm))
families <- c('Enh','EnhBiv','Het','Quies','ReprPC','TSS','TssBiv','Tx','ZNF/Rpts')
fisher_results <- lapply(c("EURv1" ,"SEX2" , "AGE"  , "BMI"), function(trait) lapply(families, function(region) my_fisher(region,trait)))
names(fisher_results) <- c("Ancestry" ,"Sex" , "Age"  , "BMI")

for (name in c("Ancestry" ,"Sex" , "Age"  , "BMI")) {
  names(fisher_results[[name]]) <- families
}

saveRDS(fisher_results, '/gpfs/projects/bsc83/Projects/GTEx_v8/Methylation/Tissues/enrichment_chromhmm_hyper_shared_CI.continous.2.rds')

# enrichment of tissue-specific positions 
sharing <- readRDS(paste0(basepath,'Projects/GTEx_v8/Methylation/Data/Sharing_DMP.rds'))
shared_cpgs <- readRDS(paste0(basepath, '/Projects/GTEx_v8/Methylation/Data/chromHMM_shared_9tissues.rds'))
families <- c('Enh','EnhBiv','Het','Quies','ReprPC','TSS','TssBiv','Tx','ZNF/Rpts')

my_fisher_nonshared <- function(type, trait, dir_value) {

  res <- sharing[sharing$trait == trait, ]
  
  # NON-SHARED set: exactly 1 tissue
  nonshared <- res[res$number == 1, ]
  
  nonshared_dir <- nonshared[nonshared$dir == dir_value, ]
  
  # universe = CpGs that appear in the sharing table for this trait (keeps ascertainment consistent)
  universe_cpgs <- unique(nonshared$CG)
  
  chrom_tissue <- shared_cpgs  # consensus chromHMM per CpG
  
  type_df<- chrom_tissue[chrom_tissue$region_chromhmm == type & chrom_tissue$name_ann %in% universe_cpgs, ]
  other_df <- chrom_tissue[chrom_tissue$region_chromhmm != type & chrom_tissue$name_ann %in% universe_cpgs, ]
  
  type_diff <- sum(type_df$name_ann %in% nonshared_dir$CG)
  type_notdiff <- nrow(type_df) - type_diff
  
  other_diff <- sum(other_df$name_ann %in% nonshared_dir$CG)
  other_notdiff <- nrow(other_df) - other_diff
  
  m <- matrix(c(type_diff, type_notdiff,
                other_diff, other_notdiff), 2, 2, byrow = TRUE)
  m[is.na(m)] <- 0
  rownames(m) <- c(type, "Other")
  colnames(m) <- c("NonShared_DMP", "Not_NonShared_DMP")
  
  f <- fisher.test(m)
  
  # return same structure you use elsewhere: list(f=..., m=type_diff)
  return(list("f" = f, "m" = type_diff))
}


# Run for Ancestry, separately by direction (AA vs EUR)
fisher_results_nonshared <- list()

fisher_results_nonshared[["AA"]] <- lapply(families, function(region)
  my_fisher_nonshared(region, trait="EURv1", dir_value=-1)
)
names(fisher_results_nonshared[["AA"]]) <- families

fisher_results_nonshared[["EA"]] <- lapply(families, function(region)
  my_fisher_nonshared(region, trait="EURv1", dir_value=1)
)
names(fisher_results_nonshared[["EA"]]) <- families

saveRDS(fisher_results_nonshared,
        paste0(basepath,'/Projects/GTEx_v8/Methylation/Tissues/enrichment_chromhmm_ancestry_nonshared_CI.rds'))

# Also report counts (for rebuttal / figure caption)
nonshared_all <- sharing[sharing$trait=="EURv1" & sharing$number==1, ]
cat("Non-shared Ancestry DMPs total:", length(unique(nonshared_all$CG)), "\n")
cat("Non-shared dir==1:", length(unique(nonshared_all$CG[nonshared_all$dir==1])), "\n")
cat("Non-shared dir==-1:", length(unique(nonshared_all$CG[nonshared_all$dir==-1])), "\n")
