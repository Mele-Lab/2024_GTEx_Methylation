#!/usr/bin/env Rscript
# @Author: Maria Sopena Rios
# @E-mail: maria.sopena@bsc.es
# @Description: Check if DMP target methylation-sensitive TFBS
# @software version: R=4.2.

library(data.table)
library(dplyr)
library(BSgenome.Hsapiens.UCSC.hg38)
library(JASPAR2022)
library(TFBSTools)


#basepath <- "/gpfs/projects/bsc83/"
basepath <- "/Users/mariasopenar/cluster/"

# downloaded methyl-sensitive TF information form Yin et al., 2017 & Gralak et al. 2024

m_sensitive_1 <- fread(paste0(basepath, "Projects/GTEx_v8/Methylation/TFBS/methylation_sensitive_TFBS_Garlak.tsv"))
m_sensitive_2 <- fread(paste0(basepath, "Projects/GTEx_v8/Methylation/TFBS/methylation_sensitive_TFBS_Yin.tsv"))
