#!/usr/bin/env Rscript
# @Author: Maria Sopena Rios
# @E-mail: maria.sopena@bsc.es
# @Description: Identification of TFBS disruption by meQTLs
# @software version: R=4.4.0

basepath <- "//home/mariasr/cluster/"
tissues <- c("Lung", "ColonTransverse" ,"Ovary", "Prostate")

results_list <- list()

for (tissue in tissues) {
  
  cat("Processing:", tissue, "\n")
  
  # ----------------------------
  # Load saved objects
  # ----------------------------
  
  dmp <- readRDS(paste0(basepath, "/Projects/GTEx_v8/Methylation/Tissues/",  tissue, "/DMP_filtered_mQTL_Promoter_Enhancer.rds"))
  all_hits <- readRDS(paste0(basepath, "/Projects/GTEx_v8/Methylation/Data/FIMO/all_hits_", tissue, ".rds"))
  snps <- readRDS(paste0(basepath, "/Projects/GTEx_v8/Methylation/Tissues/", tissue, "/SNPs_filtered_mQTL_Promoter_Enhancer.rds"))

  disrupted_snps <- unique(all_hits[status %in% c("loss","gain"), sequence_name])
  disrupted_cpgs <- unique(snps$cpg[snps$snp %in% disrupted_snps])
  
  dmp_cpgs <- unique(dmp$cpg)
  snps_with_motif <- unique(all_hits$sequence_name)
  
  overlap_cpgs <- unique(
    snps$cpg[snps$snp %in% snps_with_motif]
  )
  
  N_total <- length(dmp_cpgs)
  N_overlap <- length(overlap_cpgs)
  N_disrupt <- length(disrupted_cpgs)
  prop_overlap <- N_overlap / N_total
  prop_disrupt <- N_disrupt / N_total
  prop_disrupt_within_overlap <- ifelse(N_overlap > 0, N_disrupt / N_overlap, NA_real_)
  
  df <- data.frame(
    tissue=tissue, 
    n_total=N_total,
    n_overlap=N_overlap,
    n_disrupted=N_disrupt,
    p_overlap=prop_overlap, 
    p_disrupted_overlap=prop_disrupt_within_overlap
  )
  results_list <- rbind(results_list, df)
  
}

# plot the proportion of 

df <- results_list %>%
  mutate(
    p_disrupt_total = p_overlap * p_disrupted_overlap,
    p_overlap_not_disrupt = p_overlap - p_disrupt_total,
    p_no_overlap = 1 - p_overlap
  )

df_long <- df %>%
  select(tissue, p_no_overlap, p_overlap_not_disrupt, p_disrupt_total) %>%
  pivot_longer(-tissue, names_to="segment", values_to="proportion") %>%
  mutate(segment = recode(segment,
                          p_no_overlap = "No TFBS overlap",
                          p_overlap_not_disrupt = "TFBS overlap - not disrupted",
                          p_disrupt_total = "TFBS overlap - disrupted"))

df_long$tissue <- factor(df_long$tissue, levels = rev(df$tissue))
df_long$segment <- factor(df_long$segment, levels=rev(c("TFBS overlap - disrupted", "TFBS overlap - not disrupted","No TFBS overlap")))


ggplot(df_long, aes(x=proportion, y=tissue, alpha=segment)) +
  geom_col(width=0.7, fill="#C49122") +
  scale_x_continuous(limits=c(0,1)) +scale_alpha_manual(values=rev(c(1, 0.6, 0.2)))+
  theme_bw() +
  xlab("Proportion of cis-driven DMP CpGs") +
  ylab("") +
  theme(legend.title = element_blank(),
        axis.text.x = element_text(colour="black", size=13),
        axis.text.y = element_text(colour="black", size=14),
        legend.text = element_text(colour="black", size=13),
        axis.title.x = element_text(size=16),
        legend.spacing.y = unit(-0.05, "cm"),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.border = element_rect(colour = "black", linewidth=1))

# plot the methyl-sensitive 
library(dplyr)
library(tidyr)
library(ggplot2)

basepath <- "/Users/mariasopenar/cluster/"
tissues <- c("Lung", "ColonTransverse", "Ovary", "Prostate")

results_list <- list()

for (tissue in tissues) {
  
  prefix <- paste0(basepath,
                   "Projects/GTEx_v8/Methylation/TFBS/meth_sensitive_overlap_",
                   tissue)
  
  prop_ms  <- readRDS(paste0(prefix, "_proportion.rds"))
  dmp_cpgs <- readRDS(paste0(prefix, "_unique_cpgs_dmp.rds"))
  cpgs_ms  <- readRDS(paste0(prefix, "_cpgs_ms.rds"))
  
  N_total <- length(dmp_cpgs)
  N_overlap <- length(intersect(dmp_cpgs, cpgs_ms))
  
  df <- data.frame(
    tissue = tissue,
    p_overlap = N_overlap / N_total,
    p_no_overlap = 1 - (N_overlap / N_total)
  )
  
  results_list[[tissue]] <- df
}

df <- bind_rows(results_list)

