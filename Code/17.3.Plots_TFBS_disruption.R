#!/usr/bin/env Rscript
# @Author: Maria Sopena Rios
# @E-mail: maria.sopena@bsc.es
# @Description: Identification of TFBS disruption by meQTLs
# @software version: R=4.4.0

basepath <- "//home/mariasr/cluster/"
basepath <- "/Users/mariasopenar/cluster/"

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


p <- ggplot(df_long, aes(x=proportion, y=tissue, alpha=segment)) +
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

pdf("/Users/mariasopenar/cluster/Projects/GTEx_v8/Methylation/Plots/overlap_methyl_sensitive.pdf",
    width = 7.89, height = 3.89)
print(p)
dev.off()


# plot the methyl-sensitive 
library(dplyr)
library(tidyr)
library(ggplot2)

basepath <- "/Users/mariasopenar/cluster/"
tissues <- c("Lung", "ColonTransverse", "Ovary", "Prostate", "BreastMammaryTissue", "KidneyCortex", "Testis", "WholeBlood", "MuscleSkeletal")

results_list <- list()

for (tissue in tissues) {
  
  prefix <- paste0(basepath,
                   "Projects/GTEx_v8/Methylation/TFBS/meth_sensitive_overlap_",
                   tissue)
  
  prop_ms  <- readRDS(paste0(prefix, "proportion.rds"))
  print(prop_ms)
  dmp_cpgs <- readRDS(paste0(prefix, "_unique_cpgs_dmp.rds"))
  cpgs_ms  <- readRDS(paste0(prefix, "_cpgs_ms.rds"))
  
  N_total <- length(dmp_cpgs)
  N_overlap <- length(intersect(dmp_cpgs, cpgs_ms))
  
  df <- data.frame(
    tissue = tissue,
    p_overlap = N_overlap / N_total,
    p_no_overlap = 1 - (N_overlap / N_total),
    n_overlap=N_overlap,
    n_all=N_total
  )
  
  results_list[[tissue]] <- df
}

df <- bind_rows(results_list)

library(tidyverse)
df_long <- df[, c("p_overlap", "p_no_overlap", "tissue")] %>%
  pivot_longer(-tissue, names_to="segment", values_to="proportion") %>%
  mutate(segment = recode(segment,
                          p_overlap = "Overlap methyl-sensitive TFBS",
                          p_no_overlap = "No overlap"))

df_long$tissue <- factor(df_long$tissue, levels = rev(df$tissue))


ggplot(df_long, aes(x = proportion, y = tissue, fill = segment)) +
  geom_col(width = 0.7) +
  scale_fill_manual(values = c(
    "No overlap" = "grey",
    "Overlap methyl-sensitive TFBS" = "#777777"
  )) +
  scale_x_continuous(limits = c(0,1)) +
  theme_bw() +
  xlab("Proportion of DMP CpGs") +
  ylab("") +
  theme(
    legend.title = element_blank(),
    axis.text.x = element_text(colour="black", size=13),
    axis.text.y = element_text(colour="black", size=14),
    legend.text = element_text(colour="black", size=13),
    axis.title.x = element_text(size=16),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    panel.border = element_rect(colour = "black", linewidth = 1)
  )



# plot per trait
traits <- c("AGE","SEX2","EURv1","BMI")

df_list <- list()

for(tissue in tissues){
  for(tr in traits){
    
    prefix <- paste0(
      basepath,
      "Projects/GTEx_v8/Methylation/TFBS/meth_sensitive_overlap_",
      tissue,"_",tr
    )
    
    prop_file <- paste0(prefix,"_proportion.rds")
    dmp_file  <- paste0(prefix,"_unique_cpgs_dmp.rds")
    overlap_file <- paste0(prefix,"_cpgs_ms.rds")
    
    if(!file.exists(prop_file)) next
    
    prop <- readRDS(prop_file)
    dmp_cpgs <- readRDS(dmp_file)
    overlap <- readRDS(overlap_file)
    
    N_total <- length(dmp_cpgs)
    N_overlap <- length(overlap)
    
    df_list[[paste(tissue,tr,sep="_")]] <- data.frame(
      tissue = tissue,
      trait = tr,
      p_overlap = N_overlap/N_total,
      p_no_overlap = 1 - (N_overlap/N_total)
    )
  }
}

df <- bind_rows(df_list)
df_long <- df %>%
  pivot_longer(
    cols = c(p_overlap, p_no_overlap),
    names_to = "segment",
    values_to = "proportion"
  ) %>%
  mutate(segment = recode(segment,
                          p_overlap = "Overlap methyl-sensitive TFBS",
                          p_no_overlap = "No overlap"))

df_long$tissue <- factor(df_long$tissue, levels = rev(tissues))
df_long[df_long$trait=="AGE",]$trait <- "Age"
df_long[df_long$trait=="EURv1",]$trait <- "Ancestry"
df_long[df_long$trait=="SEX2",]$trait <- "Sex"


traits_cols <- c("Ancestry" = "#F0AE21",
                 "Age" = "#3D7CD0",
                 "Sex" =  "#3B734E",
                 "BMI" = "#CC79A7")
traits <- names(traits_cols)

p <- ggplot(df_long, aes(x = proportion, y = tissue, fill = trait, alpha=segment)) +
  geom_col(width = 0.7) +scale_alpha_manual(values=c(0.6,1))+
  facet_wrap(~trait) +
  scale_fill_manual(values=traits_cols
  ) +
  scale_x_continuous(limits=c(0,1)) +
  theme_bw() +
  xlab("Proportion of DMP CpGs") +
  ylab("") +
  theme(
    legend.title = element_blank(),
    axis.text.x = element_text(colour="black", size=13),
    axis.text.y = element_text(colour="black", size=14),
    legend.text = element_text(colour="black", size=13),
    axis.title.x = element_text(size=16),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    panel.border = element_rect(colour="black", linewidth=1)
  )

pdf("/home/mariasr/cluster/Projects/GTEx_v8/Methylation/Plots/overlap_methyl_sensitive.pdf",
    width = 11.55, height = 3.5)
print(p)
dev.off()

