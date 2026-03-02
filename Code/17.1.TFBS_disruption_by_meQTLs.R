#!/usr/bin/env Rscript
# @Author: Maria Sopena Rios
# @E-mail: maria.sopena@bsc.es
# @Description: analysis of TFBS disruption by meQTLs at DMPs.
# @software version: R=4.2.

library(argparse)
library(data.table)
library(dplyr)
library(BSgenome.Hsapiens.UCSC.hg38)
library(JASPAR2022)
library(TFBSTools)
library(seqinr)
library(data.table)
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

#ChIP Atlas annotation
#chip_atlas <- read.csv(paste0(basepath, "Projects/GTEx_v8/Methylation/TFBS/chip_seq_ALL.csv"))

# annotation positions
annotation <- read.delim(paste0(basepath, "Projects/GTEx_v8/Methylation/Data/Methylation_Epic_gene_promoter_enhancer_processed.txt"), sep = '\t', header = T)

# #read DMPs
dnp_list <-  readRDS(paste0(basepath, "/Projects/GTEx_v8/Methylation/Tissues/",tissue,"/DML_results_5_PEERs_continous.rds"))
dmp <- do.call(rbind, Map(function(df, nm) {df$trait <- nm
df$cpg <-rownames(df)
return(df)}, dnp_list, names(dnp_list)))
dmp <- dmp[dmp$adj.P.Val < 0.05,]

#mQTL data
print("Reading mQTL data")
inpath_mqtls <- "/gpfs/scratch/bsc83/MN4/bsc83/bsc83535/GTEx/v9/mQTLs/"
mqtl <- fread(paste0(inpath_mqtls, tissue, ".mQTLs.conditional.txt.gz"))

# Filter variants distance and significance
mqtl <- mqtl[mqtl$V7 < 0.05 & abs(mqtl$V3) < 250000]
mqtl$cpg <- sub(":.*", "", mqtl$V1)

#keep only cpgs that are DMPs
mqtl <- mqtl[mqtl$cpg %in% dmp$cpg]


#filter dmps
print("Filtering DMPs")
dmp <- dmp[dmp$cpg %in% mqtl$cpg, ] # that have mQTLs
dmp <- dmp[dmp$cpg %in% annotation[annotation$Type %in% c("Promoter_Associated", "Enhancer_Associated"), "IlmnID"],] # that are in enhancers and promoters

# save DMPs filtered by having an mQTL and located in promoters and enhancers
saveRDS(dmp, paste0(basepath, "/Projects/GTEx_v8/Methylation/Tissues/",tissue,"/DMP_filtered_mQTL_Promoter_Enhancer.rds"))

# Reduce set to the lead SNP per CpG
setorder(mqtl, cpg, V7)          # smallest FDR first (replace with p-value if you have it)
mqtl <- mqtl[, .SD[1], by=cpg] 

#snp parsing
print("Filtering SNPs")
snps <- unique(mqtl[mqtl$cpg %in% dmp$cpg, .(snp = V2, cpg)])
#add info
snps[, chr := sub("_.*", "", snp)]
snps[, pos := as.integer(sub("^[^_]+_([0-9]+)_.*", "\\1", snp))]
snps[, ref := sub("^[^_]+_[0-9]+_([ACGT])_.*", "\\1", snp)]
snps[, alt := sub("^[^_]+_[0-9]+_[ACGT]_([ACGT])_.*", "\\1", snp)]
snps <- snps[nchar(ref) == 1 &nchar(alt) == 1]
saveRDS(snps, paste0(basepath, "/Projects/GTEx_v8/Methylation/Tissues/",tissue,"/SNPs_filtered_mQTL_Promoter_Enhancer.rds"))


# ------------------------------------------------------------------
# Generate FASTA sequences
# ------------------------------------------------------------------

print("Generating FASTA...")
get_snp_seq <- function(chr, pos, allele, flank = 25) {
  seq <- getSeq(BSgenome.Hsapiens.UCSC.hg38,
                names = chr,
                start = pos - flank,
                end = pos + flank)
  seq <- as.character(seq)
  substr(seq, flank + 1, flank + 1) <- allele
  print(allele)
  as.character(seq)
}

param <- MulticoreParam(workers = 16)

ref_seqs <- unlist( bplapply(seq_len(nrow(snps)), function(i) get_snp_seq(snps$chr[i], snps$pos[i], snps$ref[i]),BPPARAM = param), use.names = FALSE)
alt_seqs <- unlist(bplapply(seq_len(nrow(snps)),function(i) get_snp_seq(snps$chr[i], snps$pos[i], snps$alt[i]),BPPARAM = param),use.names = FALSE)

ref_fa <- paste0("/gpfs/scratch/bsc83/MN4/bsc83/bsc83535/GTEx/v9/FIMO/", tissue, "_ref.fa")
alt_fa <- paste0("/gpfs/scratch/bsc83/MN4/bsc83/bsc83535/GTEx/v9/FIMO/", tissue, "_alt.fa")

write.fasta(as.list(ref_seqs),names=snps$snp, ref_fa, open = "w", nbchar = 60, as.string = FALSE)
write.fasta(as.list(alt_seqs),names=snps$snp, alt_fa, open = "w", nbchar = 60, as.string = FALSE)

# ------------------------------------------------------------------
# Run FIMO
# ------------------------------------------------------------------
print("Running FIMO...")
motif_file <- "/gpfs/scratch/bsc83/MN4/bsc83/bsc83535/GTEx/v9/JASPAR2022_CORE_vertebrates_non-redundant_pfms_meme.txt"
outpath_ref <- paste0("/gpfs/scratch/bsc83/MN4/bsc83/bsc83535/GTEx/v9/FIMO/", tissue, "_ref")
outpath_alt <- paste0("/gpfs/scratch/bsc83/MN4/bsc83/bsc83535/GTEx/v9/FIMO/", tissue, "_alt")

system(paste( "fimo --thresh 1e-4 --oc ",outpath_ref,  motif_file, ref_fa))
system(paste( "fimo --thresh 1e-4 --oc",outpath_alt, motif_file, alt_fa))

print("FIMO done")


# ------------------------------------------------------------------
# Process FIMO results
# ------------------------------------------------------------------

print("Processing FIMO results")
fimo_ref <- read.table(paste0("/gpfs/scratch/bsc83/MN4/bsc83/bsc83535/GTEx/v9/FIMO/", tissue, "_ref/fimo.tsv"), header = T)
fimo_alt <- read.table(paste0("/gpfs/scratch/bsc83/MN4/bsc83/bsc83535/GTEx/v9/FIMO/", tissue, "_alt/fimo.tsv"), header = T)
setDT(fimo_ref)
setDT(fimo_alt)

# Keep only hits overlapping SNP position (center of 51bp window)
center <- 26
fimo_ref <- fimo_ref[q.value < 0.05 & start <= center & stop >= center]
fimo_alt <- fimo_alt[q.value < 0.05 & start <= center & stop >= center]

# unique SNP–motif pairs
ref_hits <- unique(fimo_ref[, .(sequence_name, motif_id)])
alt_hits <- unique(fimo_alt[, .(sequence_name, motif_id)])

ref_hits[, snp_ref := 1]
alt_hits[, snp_alt := 1]

# get the motifs are lost, unchanged or gained
all_hits <- merge(
  ref_hits, alt_hits,
  by = c("sequence_name", "motif_id"),
  all = TRUE,
  suffixes = c("_ref", "_alt")
)

all_hits[, status := fifelse(
  !is.na(snp_ref) & is.na(snp_alt), "loss",
  fifelse(is.na(snp_ref) & !is.na(snp_alt), "gain", "unchanged")
)]
saveRDS(all_hits, paste0(basepath, "/Projects/GTEx_v8/Methylation/Data/FIMO/all_hits_", tissue, ".rds"))


# ------------------------------------------------------------------
# Calculate proportion of disrupted CpGs
# ------------------------------------------------------------------

disrupted_snps <- unique(all_hits[status %in% c("loss","gain"), sequence_name])
disrupted_cpgs <- unique(snps$cpg[snps$snp %in% disrupted_snps])
prop_disrupted <- length(disrupted_cpgs) / length(unique(dmp$cpg))
cat("Proportion of disrupted DMP CpGs in tissue", tissue, "=", prop_disrupted, "\n")

# save outputs
saveRDS(disrupted_cpgs, paste0(basepath, "/Projects/GTEx_v8/Methylation/Data/FIMO/disrupted_cpgs_", tissue, ".rds"))




