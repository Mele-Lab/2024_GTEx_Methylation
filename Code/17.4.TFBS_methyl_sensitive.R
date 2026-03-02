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
#library(IlluminaHumanMethylationEPICv2anno.20a1.hg38)


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
anno <- fread(paste0(scratch, "/Oliva/GPL21145_MethylationEPIC_15073387_v-1-0_processed_jose.csv"))
# adapt column names if necessary
setnames(anno,
         old = c("IlmnID","CHR","MAPINFO"),
         new = c("cpg","chr","pos"))

anno <- anno[, .(cpg, chr, pos)]

dmp <- merge(dmp, anno, by="cpg", all.x=TRUE)

cat("Missing coordinates:", sum(is.na(dmp$chr)), "\n")

dmp <- dmp[!is.na(chr) & !is.na(pos), ]

dmp_cpgs <- unique(dmp$cpg)
cat("Unique DMP CpGs:", length(dmp_cpgs), "\n")

# ----------------------------
# Define background CpGs (promoter+enhancer CpGs with coords)
# ----------------------------
bg <- data.table(cpg = unique(valid_cpgs)) # check
bg <- merge(bg, anno, by.x="cpg", by.y="IlmnID", all.x=TRUE)
bg <- bg[!is.na(chr) & !is.na(pos), ]
bg_cpgs <- unique(bg$cpg)
cat("Background CpGs (prom+enh):", length(bg_cpgs), "\n")

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


# ----------------------------
# Run FIMO
# ----------------------------
print("Running FIMO...")
motif_file <- paste0( scratch, "JASPAR2022_CORE_vertebrates_non-redundant_pfms_meme.txt")
outdir <- paste0(scratch, "FIMO/", tissue, "_meth_sensitive")
system(paste( "fimo --thresh 1e-4 --oc", outdir, motif_file, fa_file ))

fimo <- fread(paste0(outdir, "/fimo.tsv"))
center <- 26
fimo <- fimo[q.value < 0.05 & start <= center & stop >= center]
cat("Total TFBS hits:", nrow(fimo), "\n")


# ----------------------------
# Load methylation-sensitive TF lists
# ----------------------------
m_sensitive_1 <- fread(paste0(basepath, "Projects/GTEx_v8/Methylation/TFBS/methylation_sensitive_TFBS_Garlak.tsv"))
keep_classes <- c("methyl minus", "weak methyl minus", "methyl plus", "weak methyl plus")
m_sensitive_1_filt <- m_sensitive_1[classification %in% keep_classes]

m_sensitive_2 <- fread(paste0(basepath, "Projects/GTEx_v8/Methylation/TFBS/methylation_sensitive_TFBS_Yin.tsv"))
keep_classes <- c("MethylMinus and MethylPlus" , "MethylPlus", "MethylMinus", "MethylPlus and MethylMinus")
m_sensitive_2_filt <- m_sensitive_2[Call %in% keep_classes]

# get list of metht-sense-TF
meth_sensitive_tfs <- unique(c(m_sensitive_1[["TF"]], m_sensitive_2[["TF name"]]))
meth_sensitive_tfs <- meth_sensitive_tfs[!is.na(meth_sensitive_tfs)]
cat("Methylation-sensitive TFs loaded:", length(meth_sensitive_tfs), "\n")


# ----------------------------
# Map TF names -> JASPAR motif IDs
# ----------------------------
pfmList <- getMatrixSet(JASPAR2022, opts = list(collection="CORE", tax_group="vertebrates"))

motif_map <- data.frame(
  motif_id = names(pfmList),
  tf_name  = sapply(pfmList, function(x) x@name),
  stringsAsFactors = FALSE
)

# Match by TF name (case-insensitive)
motif_map$tf_name_upper <- toupper(motif_map$tf_name)
meth_sensitive_upper <- toupper(meth_sensitive_tfs)

meth_sensitive_motif_ids <- unique(motif_map$motif_id[motif_map$tf_name_upper %in% meth_sensitive_upper])
cat("Methylation-sensitive JASPAR motifs matched:", length(meth_sensitive_motif_ids), "\n")

# ----------------------------
# Filter FIMO hits to methylation-sensitive motifs
# ----------------------------
fimo_ms <- fimo[motif_id %in% meth_sensitive_motif_ids]
cat("Methylation-sensitive motif hits:", nrow(fimo_ms), "\n")

cpgs_ms <- unique(fimo_ms$sequence_name)
cat("Unique DMP CpGs overlapping methylation-sensitive motifs:", length(cpgs_ms), "\n")

prop_ms <- length(intersect(dmp_cpgs, cpgs_ms)) / length(dmp_cpgs)
out_prefix <- paste0(basepath, "Projects/GTEx_v8/Methylation/TFBS/meth_sensitive_overlap_", tissue)

saveRDS(prop_ms, paste0(out_prefix, "proportion.rds"))
saveRDS(dmp_cpgs, paste0(out_prefix, "_unique_cpgs_dmp.rds"))
saveRDS(cpgs_ms, paste0(out_prefix, "_cpgs_ms.rds"))

cat("DONE. Proportion of DMP CpGs overlapping methylation-sensitive motifs = ", prop_ms, "\n")
