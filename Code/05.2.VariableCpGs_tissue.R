
#!/usr/bin/env Rscript
# @Author: Maria Sopena Rios
# @E-mail: maria.sopena@bsc.es
# @Description: Identification of highly variable CpGs across individuals within a tissue 
# @software version: R=4.2.2

library(variancePartition)
library(BiocParallel)
library(dplyr)
first_dir <- "/gpfs/projects/bsc83/"


# Parsing
library(optparse)
parser <- OptionParser()
parser <- add_option(parser, opt_str=c("-t", "--tissue"), type="character",
                     dest="tissue",
                     help="Tissue")
parser <- add_option(parser, opt_str=c("-a", "--ancestry"), type="character",
                     dest="ancestry",
                     help="Ancestry continous / categorical")
parser <- add_option(parser, opt_str=c("-r", "--remove_admixed"), type="logical",
                     dest="remove_admixed",
                     help="Remove")
parser <- add_option(parser, opt_str=c("-s", "--smoking_status"), type="logical",
                     dest="smoking_status",
                     help="smoking information")
options=parse_args(parser)
tissue=options$tissue


tissues <- c("BreastMammaryTissue","ColonTransverse","KidneyCortex","Lung",
             "MuscleSkeletal","Ovary","Prostate","Testis","WholeBlood")
# files <- list.files('/gpfs/projects/bsc83/Projects/GTEx_v8/Methylation/Data/EpiMap/', pattern='.bed.gz',full.names=T)
# 
# annotation <- read.csv(paste0(first_dir, "/Projects/GTEx_v8/Methylation/Data/GPL21145_MethylationEPIC_15073387_v-1-0_processed.csv"))
# 
# names_chrom <- c("Lung", "ColonTransverse", "Ovary", "Prostate", "Breast", "MuscleSkeletal", "KidneyCortex", "Testis", "PBMC")
# chromhmm <- lapply(names_chrom, function(tis)
#   read.delim(files[grep(tis, files)], sep='\t', header=F))
# names(chromhmm) <- tissues
# 
# ann_bed <- annotation[!is.na(annotation$MAPINFO) & !is.na(annotation$CHR),] %>%
#   dplyr::select(chrom=CHR, start=MAPINFO, end=MAPINFO, name=IlmnID) %>% distinct()
# ann_bed$chrom <- paste0('chr',ann_bed$chrom)
# head(ann_bed)
# ann_bed$start <- ann_bed$start-1
# 
# library(valr)
# chromhmm_cpgs <- lapply(tissues, function(tis) {
#   chrom_df <- chromhmm[[tis]][,c(1:4)]
#   colnames(chrom_df) <- c('chrom','start','end','region')
#   bed_intersect(ann_bed, chrom_df, suffix = c("_ann", "_chromhmm"))})
# names(chromhmm_cpgs) <- tissues
# 
# print(tissue)
# meta <- readRDS(paste0(first_dir, "/Projects/GTEx_v8/Methylation/Tissues/", tissue, "/metadata.rds"))
# data <- readRDS(paste0(first_dir, "/Projects/GTEx_v8/Methylation/Tissues/", tissue, "/data.rds"))
# beta <- data
# # align
# rownames(meta) <- meta$SUBJID
# common <- intersect(colnames(beta), rownames(meta))
# beta <- beta[, common, drop=FALSE]
# beta <- sapply(beta, as.numeric)
# rownames(beta) <- rownames(data)
# meta <- meta[common, , drop=FALSE]
# print("Getting highly variable CpG...")
# # get residuals 
# resid_beta <- limma::removeBatchEffect(
#   beta,
#   covariates = as.matrix(meta[, c("PEER1","PEER2","PEER3","PEER4","PEER5")])
# )
# 
# # get variability of beta values 
# var_cpg <- apply(resid_beta, 1, var, na.rm=TRUE)
# 
# # define highly variable CpGs as being at the top 5% 
# thr <- quantile(var_cpg, 0.95, na.rm=TRUE)
# high_var_cpgs <- names(var_cpg)[var_cpg >= thr]
# #saveRDS(high_var_cpgs, paste0("/gpfs/projects/bsc83/Projects/GTEx_v8/Methylation/varPart/", tissue, "_highly_variable_CpGs.rds"))
# saveRDS(names(var_cpg), paste0("/gpfs/projects/bsc83/Projects/GTEx_v8/Methylation/varPart/", tissue, "_all_CpGs.rds"))


# get the enrichments in Enhancers and Promoters 


anno <- read.delim(paste0("/gpfs/projects/bsc83/Projects/GTEx_v8/Methylation/Data/Methylation_Epic_gene_promoter_enhancer_processed.txt"), sep = '\t', header = T)

high_var_cpgs <- readRDS(paste0("/gpfs/projects/bsc83/Projects/GTEx_v8/Methylation/varPart/", tissue, "_highly_variable_CpGs.rds"))
var_cpg <-readRDS(paste0("/gpfs/projects/bsc83/Projects/GTEx_v8/Methylation/varPart/", tissue, "_all_CpGs.rds"))

my_fisher_broad <- function(type, tissue, variable_cpg, universe){
  anno_universe <- intersect(universe, unique(anno$IlmnID))
  var_set <- intersect(variable_cpgs, anno_universe)
  # CpGs in the given state
  in_type <- intersect(anno_universe, anno$IlmnID[anno$Type == type])
  
  # 2x2 counts
  a <- sum(var_set %in% in_type)                 # variable & in_type
  b <- length(in_type) - a                       # not variable & in_type
  c <- length(var_set) - a                       # variable & not in_type
  d <- (length(anno_universe) - length(in_type)) - c   # not variable & not in_type

  ### test significance
  m <- matrix(c(a,b,c,d), nrow=2, byrow=TRUE)
  rownames(m) <- c("InState", "NotInState")
  colnames(m) <- c("Variable", "NotVariable")
  print(m)

  m[is.na(m)] <- 0
  #m <- m[c(type,paste0('No ',type)),]
  rownames(m) <- c(type, "Other")
  colnames(m) <- c("Variable","No_Variable")
  print(m)
  f <- fisher.test(m)
  print(f)
  return(list("f" = f, "m" = m))
  
}

types <- c("Gene_Associated", "Enhancer_Associated", "Promoter_Associated")
fisher_results <- lapply(types, function(type) my_fisher_broad(type,tissue,high_var_cpgs, var_cpg ))
names(fisher_results) <-types

saveRDS(fisher_results, paste0("/gpfs/projects/bsc83/Projects/GTEx_v8/Methylation/varPart/", tissue, "_highly_variable_CpGs_enrichment_BroadClassification.rds"))




# 
# my_fisher <- function(type, tissue, variable_cpgs, universe){
#   
#   chrom_tissue <- chromhmm_cpgs[[tissue]]
#   chrom_tissue$region_chromhmm_new <- chrom_tissue$region_chromhmm
#   chrom_tissue <- chrom_tissue[chrom_tissue$.overlap ==1,]
#   # 
#   chrom_tissue$region_chromhmm_new[chrom_tissue$region_chromhmm %in% c("TssFlnkD", "TssFlnk", "TssFlnkU","TssA")] <- "TSS"
#   chrom_tissue$region_chromhmm_new[chrom_tissue$region_chromhmm %in% c("EnhA2", "EnhA1","EnhWk","EnhG1", "EnhG2")] <- "Enh"
#   chrom_tissue$region_chromhmm_new[chrom_tissue$region_chromhmm %in% c("ReprPCWk","ReprPC")] <- "ReprPC"
#   chrom_tissue$region_chromhmm_new[chrom_tissue$region_chromhmm %in% c("TxWk","Tx")] <- "Tx"
#   
#   # Universe must be CpGs that have a chromHMM assignment in this tissue
#   anno_universe <- intersect(universe, unique(chrom_tissue$name_ann))
#   var_set <- intersect(variable_cpgs, anno_universe)
#   
#   # CpGs in the given state
#   in_type <- intersect(anno_universe, chrom_tissue$name_ann[chrom_tissue$region_chromhmm_new == type])
#  
#   # CpGs in the given state
#   in_type <- intersect(anno_universe, chrom_tissue$name_ann[chrom_tissue$region_chromhmm_new == type])
#   
#   # 2x2 counts
#   a <- sum(var_set %in% in_type)                 # variable & in_type
#   b <- length(in_type) - a                       # not variable & in_type
#   c <- length(var_set) - a                       # variable & not in_type
#   d <- (length(anno_universe) - length(in_type)) - c   # not variable & not in_type
#   
#   ### test significance
#   m <- matrix(c(a,b,c,d), nrow=2, byrow=TRUE)
#   rownames(m) <- c("InState", "NotInState")
#   colnames(m) <- c("Variable", "NotVariable")
#   print(m)
#   
#   m[is.na(m)] <- 0
#   #m <- m[c(type,paste0('No ',type)),]
#   rownames(m) <- c(type, "Other")
#   colnames(m) <- c("Variable","No_Variable")
#   print(m)
#   f <- fisher.test(m)
#   print(f)
#   return(list("f" = f, "m" = m))
#   
# }
# 
# print("Running fisher...")
# # Two-tailed Fisher test
# #families <- as.vector(unique(shared_cpgs$region_chromhmm))
# families <- c('Enh','EnhBiv','Het','Quies','ReprPC','TSS','TssBiv','Tx','ZNF/Rpts')
# fisher_results <- lapply(families, function(type) my_fisher(type,tissue,high_var_cpgs, rownames(data) ))
# names(fisher_results) <-families
# 
# saveRDS(fisher_results, paste0("/gpfs/projects/bsc83/Projects/GTEx_v8/Methylation/varPart/", tissue, "_highly_variable_CpGs_enrichment_chromHMM.rds"))
# 

# get the functional enrichments

library(missMethyl)

GOenrichments <- list()
for (tissue in tissues) {
  print(tissue)
  GOenrichments[[tissue]] <- list()
  high_var_cpgs <- readRDS(paste0(first_dir, "/Projects/GTEx_v8/Methylation/varPart/", tissue, "_highly_variable_CpGs.rds"))
  all_cpgs <- readRDS(paste0(first_dir, "/Projects/GTEx_v8/Methylation/varPart/", tissue, "_all_CpGs.rds"))
  res <- missMethyl::gometh(high_var_cpgs, all.cpg=all_cpgs,
                            collection="GO", array.type="EPIC")
  res <- res[res$ONTOLOGY=="BP",]
  print(table(res$FDR<0.05))
  saveRDS(res,paste0("/gpfs/projects/bsc83/Projects/GTEx_v8/Methylation/varPart/", tissue, "_highly_variable_CpGs_Functiona_enrichment.rds") )
  
  if (sum(res$FDR<0.05) > 0) {
    GOenrichments[[tissue]] <- res[res$FDR<0.05,]
  }
}

saveRDS(GOenrichments,paste0("/gpfs/projects/bsc83/Projects/GTEx_v8/Methylation/varPart/All_tissues_highly_variable_CpGs_Functiona_enrichment.rds") )



