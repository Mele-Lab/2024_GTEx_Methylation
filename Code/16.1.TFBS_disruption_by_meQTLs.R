#!/usr/bin/env Rscript
# @Author: Maria Sopena Rios
# @E-mail: maria.sopena@bsc.es
# @Description: analysis of TFBS disruption by meQTLs at DMPs.
# @software version: R=4.2.

library(argparse)
library(data.table)
library(dplyr)
library(BSgenome.Hsapiens.UCSC.hg38)
basepath <- "/gpfs/projects/bsc83/"

shhh <- suppressPackageStartupMessages
shhh(library(optparse))
option_list = list(
  make_option(c("--tissue"), action="store", default=NA, type='character'))
opt = parse_args(OptionParser(option_list=option_list))

tissue <- opt$tissue 

cat("Processing tissue:", tissue, "\n")
cat(Sys.time(), "\n")
#tissue <- "Lung"

#ChIP Atlas annotation
#chip_atlas <- read.csv(paste0(basepath, "Projects/GTEx_v8/Methylation/TFBS/chip_seq_ALL.csv"))

# annotation positions
annotation <- read.delim(paste0(basepath, "Projects/GTEx_v8/Methylation/Data/Methylation_Epic_gene_promoter_enhancer_processed.txt"), sep = '\t', header = T)

#read DMPs
dnp_list <-  readRDS(paste0(basepath, "/Projects/GTEx_v8/Methylation/Tissues/",tissue,"/DML_results_5_PEERs_continous.rds"))
dmp <- do.call(rbind, Map(function(df, nm) {df$trait <- nm
df$cpg <-rownames(df)
return(df)}, dnp_list, names(dnp_list)))


#mQTL data
print("Reading mQTL data")
inpath_mqtls <- "/gpfs/scratch/bsc83/MN4/bsc83/bsc83535/GTEx/v9/mQTLs/"
mqtl <- fread(paste0(inpath_mqtls, tissue, ".mQTLs.conditional.txt.gz"))
# Filter variants distance and significance
mqtl <- mqtl[mqtl$V7 < 0.05 & abs(mqtl$V3) < 250000]
mqtl$cpg <- sub(":.*", "", mqtl$V1)


#filter dmps
dmp <- dmp[dmp$cpg %in% mqtl$cpg, ] # that have mQTLs
dmp <- dmp[dmp$cpg %in% annotation[annotation$Type %in% c("Promoter_Associated", "Enhancer_Associated"), "IlmnID"],] # that are in enhancers and promoters

#snp parsing
snps <- unique(mqtl[mqtl$cpg %in% dmp$cpg, .(snp = V2, cpg)])

snps[, chr := sub(":.*", "", snp)]
snps[, pos := as.integer(sub(".*:(\\d+):.*", "\\1", snp))]
snps[, ref := sub(".*:(.):.*", "\\1", snp)]
snps[, alt := sub(".*:(.).*", "\\1", snp)]





