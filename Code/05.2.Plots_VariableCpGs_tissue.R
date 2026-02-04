
#!/usr/bin/env Rscript
# @Author: Maria Sopena Rios
# @E-mail: maria.sopena@bsc.es
# @Description: Plot enrichemnt of highly variable CpGs across individuals within a tissue 
# @software version: R=4.2.2


first_dir <- "/Users/mariasopenar/cluster/"

tissues <- c("BreastMammaryTissue", "ColonTransverse" ,"KidneyCortex", "Lung", "MuscleSkeletal" ,"Ovary", "Prostate", "Testis", "WholeBlood")

fisher_results <- lapply(tissues, function(tissue) 
  readRDS(paste0(first_dir, "/Projects/GTEx_v8/Methylation/varPart/", tissue, "_highly_variable_CpGs_enrichment_chromHMM.rds")))
