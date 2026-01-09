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
library(seqinr)
library(IlluminaHumanMethylationEPICv2anno.20a1.hg38)


#basepath <- "/gpfs/projects/bsc83/"
basepath <- "/gpfs/projects/bsc83/"
scratch  <- "/gpfs/scratch/bsc83/MN4/bsc83/bsc83535/GTEx/v9/"


shhh <- suppressPackageStartupMessages
shhh(library(optparse))
option_list = list(
  make_option(c("--tissue"), action="store", default=NA, type='character'))
opt = parse_args(OptionParser(option_list=option_list))

tissue <- opt$tissue 
cat("Processing tissue:", tissue, "\n")
cat(Sys.time(), "\n")

# annotation positions
annotation <- read.delim(paste0(basepath, "Projects/GTEx_v8/Methylation/Data/Methylation_Epic_gene_promoter_enhancer_processed.txt"), sep = '\t', header = T)
anno <- read.csv(past0(scratch, "/Oliva/GPL21145_MethylationEPIC_15073387_v-1-0_processed_jose.csv"))
#read DMPs
dnp_list <-  readRDS(paste0(basepath, "/Projects/GTEx_v8/Methylation/Tissues/",tissue,"/DML_results_5_PEERs_continous.rds"))
dmp <- do.call(rbind, Map(function(df, nm) {df$trait <- nm
df$cpg <-rownames(df)
return(df)}, dnp_list, names(dnp_list)))
dmp <- dmp[dmp$adj.P.Val < 0.05,]

# keep enhancers and promoters
dmp <- dmp[dmp$cpg %in% annotation[annotation$Type %in% c("Promoter_Associated", "Enhancer_Associated"), "IlmnID"],] # that are in enhancers and promoters
cat("Total DMPs:", nrow(dmp), "\n")


# add chromosomal location in DMP results 
anno <- getAnnotation(IlluminaHumanMethylationEPICv2anno.20a1.hg38)
anno <- as.data.table(anno, keep.rownames = "IlmnID")
dmp <- merge(dmp, anno[, .(IlmnID, chr, pos)], by = "IlmnID", all.x = TRUE)

cat("Missing coordinates:", sum(is.na(dmp$chr)), "\n")

#generate FASTA sequences to input to FIMO 
get_snp_seq <- function(chr, pos, allele, flank = 25) {
  seq <- getSeq(BSgenome.Hsapiens.UCSC.hg38,
                names = chr,
                start = pos - flank,
                end = pos + flank)
  as.character(seq)
  return(seq)
}

fa_file<- paste0("/gpfs/scratch/bsc83/MN4/bsc83/bsc83535/GTEx/v9/FIMO/", tissue, ".fa")

seqs <- mapply(get_snp_seq, dmp$chr, dmp$pos)
write.fasta( sequences = as.list(seqs),  names = dmp$cpg,  file.out = fa_file, nbchar = 60 )

motif_file <- paste0( scratch, "JASPAR2022_CORE_vertebrates_non-redundant_pfms_meme.txt"
)

outdir <- paste0(scratch, "FIMO/", tissue, "_meth_sensitive")

system(paste( "fimo --thresh 1e-4 --oc", outdir, motif_file, fa_file ))


fimo <- fread(paste0(outdir, "/fimo.tsv"))
fimo <- fimo[q.value < 0.05]

cat("Total TFBS hits:", nrow(fimo), "\n")


# downloaded methyl-sensitive TF information form Yin et al., 2017 & Gralak et al. 2024

m_sensitive_1 <- fread(paste0(basepath, "Projects/GTEx_v8/Methylation/TFBS/methylation_sensitive_TFBS_Garlak.tsv"))
m_sensitive_2 <- fread(paste0(basepath, "Projects/GTEx_v8/Methylation/TFBS/methylation_sensitive_TFBS_Yin.tsv"))
