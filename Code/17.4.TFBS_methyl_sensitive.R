#!/usr/bin/env Rscript
# @Author: Maria Sopena Rios
# @E-mail: maria.sopena@bsc.es
# @Description: Check if DMP target methylation-sensitive TFBS
# @software version: R=4.2.

library(data.table)
library(dplyr)
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


# # ---------------------------------------------------------
# # Load annotation (same one used in your working script)
# # ---------------------------------------------------------
# 
annotation <- fread(paste0(basepath,
         "Projects/GTEx_v8/Methylation/Data/Methylation_Epic_gene_promoter_enhancer_processed.txt"))

annotation <- annotation[annotation$Type %in% c("Promoter_Associated","Enhancer_Associated"),]


# extract chr and pos from Phantom5_Enhancers column
annotation[, chr := sub(":.*", "", Phantom5_Enhancers)]
annotation[, pos := as.integer(sub(".*:(\\d+)-.*", "\\1", Phantom5_Enhancers))]

#read DMPs
dnp_list <-  readRDS(paste0(basepath, "/Projects/GTEx_v8/Methylation/Tissues/",tissue,"/DML_results_5_PEERs_continous.rds"))
dmp <- do.call(rbind, Map(function(df, nm) {df$trait <- nm
df$cpg <-rownames(df)
return(df)}, dnp_list, names(dnp_list)))
dmp <- dmp[dmp$adj.P.Val < 0.05,]

# keep enhancers and promoters
dmp <- dmp[dmp$cpg %in% annotation$IlmnID, ]
cat("Total DMP rows:", nrow(dmp), "\n")

# ---------------------------------------------------------
# Add CpG genomic coordinates
# ---------------------------------------------------------

dmp <- merge(
  dmp,
  annotation[, .(IlmnID, chr, pos)],
  by.x="cpg",
  by.y="IlmnID"
)

dmp <- as.data.table(dmp)
# ensure chr format matches BSgenome
dmp[, chr := paste0("chr", gsub("chr","",chr))]

valid_chr <- paste0("chr", c(1:22,"X","Y"))
dmp <- dmp[dmp$chr %in% valid_chr,]
# 
# # ----------------------------
# # Define background CpGs (promoter+enhancer CpGs with coords)
# # ----------------------------
# # ---------------------------------------------------------
# # Generate FASTA around CpG
# # ---------------------------------------------------------
# 
# cat("Generating FASTA...\n")
# 
# library(BSgenome.Hsapiens.UCSC.hg38)
# genome <- BSgenome.Hsapiens.UCSC.hg38
# chr_lengths <- seqlengths(genome)
# 
# get_cpg_seq <- function(chr,pos,flank=25){
#   
#   if(!(chr %in% names(chr_lengths))) return(NA_character_)
#   
#   start <- max(1, pos - flank)
#   end <- min(chr_lengths[[chr]], pos + flank)
#   
#   if(start > end) return(NA_character_)
#   
#   as.character(getSeq(genome, chr, start=start, end=end))
# }
# 
# seqs <- vapply(
#   seq_len(nrow(dmp)),
#   function(i) get_cpg_seq(dmp$chr[i], dmp$pos[i]),
#   character(1)
# )
# 
# valid <- !is.na(seqs)
# 
# dmp <- dmp[valid]
# seqs <- seqs[valid]
# 
# fa_file <- paste0(scratch,"FIMO/",tissue,"_CpG.fa")
# write.fasta( sequences = as.list(seqs),  names = dmp$cpg,  file.out = fa_file, nbchar = 60 )
# 
# # ----------------------------
# # Run FIMO
# # ----------------------------
# print("Running FIMO...")
# motif_file <- paste0( scratch, "JASPAR2022_CORE_vertebrates_non-redundant_pfms_meme.txt")
outdir <- paste0(scratch, "FIMO/", tissue, "_meth_sensitive")
# system(paste( "fimo --thresh 1e-4 --oc", outdir, motif_file, fa_file ))

fimo <- fread(paste0(outdir, "/fimo.tsv"))
center <- 26
fimo <- fimo[fimo$`q-value` < 0.05 & fimo$start <= center & fimo$stop >= center]
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

traits <- unique(dmp$trait)

for(tr in traits){
  
  cat("Processing trait:", tr, "\n")
  
  dmp_trait <- dmp[dmp$trait == tr, ]
  
  dmp_cpgs <- unique(dmp_trait$cpg)
  cpgs_ms <- unique(fimo_ms$sequence_name)
  
  overlap <- intersect(dmp_cpgs, cpgs_ms)
  
  prop_ms <- length(overlap) / length(dmp_cpgs)
  
  out_prefix <- paste0(
    basepath,
    "Projects/GTEx_v8/Methylation/TFBS/meth_sensitive_overlap_",
    tissue,
    "_",
    tr
  )
  
  saveRDS(prop_ms, paste0(out_prefix, "_proportion.rds"))
  saveRDS(dmp_cpgs, paste0(out_prefix, "_unique_cpgs_dmp.rds"))
  saveRDS(overlap, paste0(out_prefix, "_cpgs_ms.rds"))
  
  cat(
    "Trait:", tr,
    "- proportion overlapping methyl-sensitive TFBS =",
    prop_ms, "\n"
  )
}