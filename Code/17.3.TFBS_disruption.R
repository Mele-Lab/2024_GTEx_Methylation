#!/usr/bin/env Rscript
# @Author: Maria Sopena Rios
# @E-mail: maria.sopena@bsc.es
# @Description: Identification of TFBS disruption by meQTLs
# @software version: R=4.4.0

setwd(paste0(basepath, "/Projects/GTEx_v8/Methylation/"))

tissues <- list.dirs("Tissues/", full.names = F)[-1]
tissues <- tissues[-grep('Old',tissues)]
names_table <- cbind(tissues, "short"=c("Lng", "Dig", "ALL", "Prs", "Brs", "Kid", "ALL", "Bld", "Myo"))

for (tissue in tissues){
  annotation <- read.csv(paste0(data_path, "chip_seq_", short,".csv"))[,-1]
  colnames(annotation)[5] <- c("TF")
  print(paste("We will be testing", length(unique(annotation$TF)), "TFs"))
} 
  
results_DML <- lapply(tissues, function(tis) 
  readRDS(paste0(basepath, "/Projects/GTEx_v8/Methylation/Tissues/",tis,"/DML_results_5_PEERs_continous.rds")))
names(results_DML) <- tissues
