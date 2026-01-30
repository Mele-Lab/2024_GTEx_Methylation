#!/usr/bin/env Rscript
# @Author: Maria Sopena Rios
# @E-mail: maria.sopena@bsc.es
# @Description: Are the effects of new mutations on DMPs biased towars lower or higher methylation?
# @software version: R=4.2.

library(argparse)
library(data.table)
library(dplyr)
library(BiocParallel)


basepath <- "/gpfs/projects/bsc83/"

shhh <- suppressPackageStartupMessages
shhh(library(optparse))
option_list = list(
  make_option(c("--tissue"), action="store", default=NA, type='character'))
opt = parse_args(OptionParser(option_list=option_list))

tissue <- opt$tissue 
#tissue <- "MuscleSkeletal"

cat("Processing tissue:", tissue, "\n")
cat(Sys.time(), "\n")

# annotation cpg positions
print("Reading CpG annotation position")
cpg_anno <- read.delim(paste0(basepath, "Projects/GTEx_v8/Methylation/Data/Methylation_Epic_gene_promoter_enhancer_processed.txt"), sep = '\t', header = T)

print("Reading VCF information")
#ancestral allele annotation 
#add VCF 

#mQTL data
print("Reading mQTL data")
inpath_mqtls <- "/gpfs/scratch/bsc83/MN4/bsc83/bsc83535/GTEx/v9/mQTLs/"
mqtl <- fread(paste0(inpath_mqtls, tissue, ".mQTLs.conditional.txt.gz"))

# Filter variants distance and significance 
mqtl <- mqtl[mqtl$V7 < 0.05 & abs(mqtl$V3) < 250000]
mqtl$cpg <- sub(":.*", "", mqtl$V1)

# Reduce set to the lead SNP per CpG
setorder(mqtl, cpg, V7) 
mqtl <- mqtl[, .SD[1], by=cpg] 

#merging 
# Merge
dt <- merge(mqtl, anno, by = "snp_id", all.x = TRUE)
dt <- merge(dt, cpg_anno, by = "cpg", all.x = TRUE)

# Keep biallelic SNPs with known AA
dt <- dt[!is.na(ancestral_allele)]
dt <- dt[nchar(effect_allele) == 1 & nchar(other_allele) == 1 & nchar(ancestral_allele) == 1]


# -------------------------
# 1) Determine derived allele
# -------------------------
# derived = allele that is not ancestral, among {EA, OA}
dt[, derived_allele := fifelse(effect_allele != ancestral_allele & other_allele == ancestral_allele, effect_allele,
                               fifelse(other_allele != ancestral_allele & effect_allele == ancestral_allele, other_allele,
                                       NA_character_))]

# Drop ambiguous cases (AA doesn't match either allele; or strand issues)
dt <- dt[!is.na(derived_allele)]

# -------------------------
# 2) Compute effect of DERIVED allele on methylation
# -------------------------
# beta is defined as effect of effect_allele (EA) on methylation.
# If derived == EA: derived_beta = beta
# If derived == OA: derived_beta = -beta   (because swapping alleles flips sign)
dt[, derived_beta := fifelse(derived_allele == effect_allele, beta, -beta)]

# Direction label
dt[, derived_direction := fifelse(derived_beta > 0, "Derived increases methylation",
                                  fifelse(derived_beta < 0, "Derived decreases methylation", "Zero"))]

# Typically drop zeros
dt <- dt[derived_direction != "Zero"]

# -------------------------
# 3) Summarize per region + binomial test (bias from 50:50)
# -------------------------
summ <- dt[, .(
  n = .N,
  n_up = sum(derived_beta > 0),
  n_down = sum(derived_beta < 0),
  prop_up = mean(derived_beta > 0)
), by = region_class]

# Exact binomial test + CI per region
summ[, c("p_value", "ci_low", "ci_high") := {
  bt <- binom.test(n_up, n, p = 0.5)
  list(bt$p.value, bt$conf.int[1], bt$conf.int[2])
}, by = region_class]

# Multiple-testing correction across regions (optional)
summ[, p_adj := p.adjust(p_value, method = "BH")]

print(summ[order(p_adj)])

