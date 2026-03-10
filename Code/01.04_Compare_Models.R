# compare models 

first_dir <- "~/"
#setwd(paste0(first_dir, "marenostrum/Projects/GTEx_v8/Methylation/"))
basepath <- "/Users/mariasopenar/cluster/"
project_path <- paste0(basepath, "/Projects/GTEx_v8/Methylation/")

tissues <- c("Lung", "ColonTransverse", "Ovary", "Prostate", "BreastMammaryTissue", "MuscleSkeletal", "KidneyCortex", "Testis", "WholeBlood")

DML_main <- lapply(tissues, function(tissue) readRDS(paste0(project_path, "/Tissues/", tissue, "/DML_results_5_PEERs_continous.rds")))
names(DML_main) <- tissues

DML_new<- lapply(tissues, function(tissue) readRDS(paste0(project_path, "/Tissues/", tissue, "/DML_results_5_PEERs_Ancestry_continous_remove_admixed_FALSE_smoking_FALSE.peer.rds")))


DML_admixed<- lapply(tissues, function(tissue) readRDS(paste0(project_path, "/Tissues/", tissue, "/DML_results_5_PEERs_Ancestry_continous_remove_admixed_TRUE_smoking_FALSE.peer.rds")))
names(DML_admixed) <- tissues
DML_smoking<- lapply(tissues, function(tissue) readRDS(paste0(project_path, "/Tissues/", tissue, "/DML_results_5_PEERs_Ancestry_continous_remove_admixed_FALSE_smoking_TRUE.peer.rds")))
names(DML_smoking) <- tissues
DML_categorical<- lapply(tissues, function(tissue) readRDS(paste0(project_path, "/Tissues/", tissue, "/DML_results_5_PEERs_Ancestry_categorical_remove_admixed_FALSE_smoking_FALSE.peer.rds")))
names(DML_categorical) <- tissues
for (tis in tissues) {
  
  cat("Processing:", tis, "\n")
  
  # Full analysis
  res_full_tis <- DML_smoking[[tis]][["EURv1"]]
  
  all_full <- rownames(res_full_tis)
  dmp_full <- rownames(res_full_tis[res_full_tis$adj.P.Val < 0.05, ])
  
  saveRDS(all_full,
          paste0(project_path, "/Tissues/", tis,"/Ancestry_continous_smoking_tested.rds" ))
  
  saveRDS(dmp_full,
          paste0(project_path, "/Tissues/", tis,"/Ancestry_continous_smoking_dmp.rds"))
}

for (tis in tissues) {
  
  cat("Processing:", tis, "\n")
  
  # Full analysis
  res_full_tis <- DML_admixed[[tis]][["EURv1"]]
  
  all_full <- rownames(res_full_tis)
  dmp_full <- rownames(res_full_tis[res_full_tis$adj.P.Val < 0.05, ])
  
  saveRDS(all_full,
          paste0(project_path, "/Tissues/", tis,"/Ancestry_continous_admixed_tested.rds" ))
  
  saveRDS(dmp_full,
          paste0(project_path, "/Tissues/", tis,"/Ancestry_continous_admixed_dmp.rds"))
}


for (tis in tissues) {
  
  cat("Processing:", tis, "\n")
  
  # Full analysis
  res_full_tis <- DML_main[[tis]][["EURv1"]]
  
  all_full <- rownames(res_full_tis)
  dmp_full <- rownames(res_full_tis[res_full_tis$adj.P.Val < 0.05, ])
  
  saveRDS(all_full,
          paste0(project_path, "/Tissues/", tis,"/Ancestry_continous_tested.rds" ))
  
  saveRDS(dmp_full,
          paste0(project_path, "/Tissues/", tis,"/Ancestry_continous_dmp.rds"))
}



for (tis in tissues) {
  
  cat("Processing:", tis, "\n")
  
  # Full analysis
  res_full_tis <- DML_categorical[[tis]][["EURv1"]]
  
  all_full <- rownames(res_full_tis)
  dmp_full <- rownames(res_full_tis[res_full_tis$adj.P.Val < 0.05, ])
  
  saveRDS(all_full,
          paste0(project_path, "/Tissues/", tis,"/Ancestry_categorical_tested.rds" ))
  
  saveRDS(dmp_full,
          paste0(project_path, "/Tissues/", tis,"/Ancestry_categorical_dmp.rds"))
}

# perform fisher enrichments 

compute_overlap_light <- function(tis, path_all,path_dmp) {
  
  all_full <- readRDS( paste0(project_path, "/Tissues/", tis,"/Ancestry_continous_tested.rds" ))
  dmp_full <- readRDS( paste0(project_path, "/Tissues/", tis,"/Ancestry_continous_dmp.rds" ))
  
  all_noadmix <- readRDS(path_all)
  dmp_noadmix <- readRDS(path_dmp)
  
  universe <- intersect(all_full, all_noadmix)
  
  dmp_full <- intersect(dmp_full, universe)
  dmp_noadmix <- intersect(dmp_noadmix, universe)
  
  overlap <- length(intersect(dmp_full, dmp_noadmix))
  only_full <- length(setdiff(dmp_full, dmp_noadmix))
  only_noadmix <- length(setdiff(dmp_noadmix, dmp_full))
  neither <- length(universe) - overlap - only_full - only_noadmix
  
  contingency <- matrix(c(overlap,
                          only_full,
                          only_noadmix,
                          neither),
                        nrow = 2)
  
  fisher_res <- fisher.test(contingency)
  
  data.frame(
    tissue = tis,
    oddsRatio = as.numeric(fisher_res$estimate),
    CI_down = fisher_res$conf.int[1],
    CI_up = fisher_res$conf.int[2],
    pvalue = fisher_res$p.value,
    sig = fisher_res$p.value < 0.05
  )
}

results_overlap_adimixed <- do.call(rbind,
                           lapply(tissues, function(tis) compute_overlap_light(tis, paste0(project_path, "/Tissues/", tis,"/Ancestry_continous_admixed_tested.rds"), paste0(project_path, "/Tissues/", tis,"/Ancestry_continous_admixed_dmp.rds"))))



results_overlap_categorical <- do.call(rbind,
                                    lapply(tissues, function(tis) compute_overlap_light(tis, paste0(project_path, "/Tissues/", tis,"/Ancestry_categorical_tested.rds"), paste0(project_path, "/Tissues/", tis,"/Ancestry_categorical_dmp.rds"))))


results_overlap_smoking <- do.call(rbind,
                                       lapply(tissues, function(tis) compute_overlap_light(tis, paste0(project_path, "/Tissues/", tis,"/Ancestry_continous_smoking_tested.rds"), paste0(project_path, "/Tissues/", tis,"/Ancestry_continous_smoking_dmp.rds"))))



#fisher_results <- results_overlap_adimixed
plot_overlap<- function(fisher_results){
  fisher_results$p.adjust <- p.adjust(fisher_results$pvalue, method = "BH")
  fisher_results$sig[fisher_results$p.adjust< 0.05] <- "FDR < 0.05"
  fisher_results$sig[fisher_results$p.adjust>= 0.05] <- "FDR >= 0.05"
  fisher_results$sig <- factor(fisher_results$sig, levels=c("FDR >= 0.05", "FDR < 0.05"))
  fisher_results$tissue <- factor(fisher_results$tissue, levels=rev(tissues))

  g1 <- ggplot(fisher_results, aes(x=log2(oddsRatio), y=tissue, alpha=sig)) +
    geom_errorbar(aes(xmin=log2(CI_down), xmax=log2(CI_up)), width=.3,  color="#87409c") +
    geom_vline(xintercept = 0) +
    #xlim(0,20) + #Only for Lung to show the 0
    geom_point(size=3, color="#87409c") + ylab('') + theme_bw() +
    #scale_colour_manual(values=colors_traits[[trait]]) +
    xlab("log2(Odds ratio)") +
    scale_alpha_discrete(range = c(0.4, 1), drop = FALSE) +
    theme(legend.title = element_blank(),
          axis.text.x = element_text(colour="black", size=12),
          axis.text.y = element_text(colour="black", size=12),
          legend.text = element_text(colour="black", size=12),
          axis.title.x = element_text(size=12),
          legend.spacing.y = unit(-0.05, "cm"),
          panel.grid.major = element_blank(),
          panel.grid.minor = element_blank(),
          panel.border = element_rect(colour = "black", linewidth=1), legend.position="top") +
    scale_y_discrete(breaks=tissues)# + xlim(0, 3)
  return(g1)
}

p_admixed <- plot_overlap(results_overlap_adimixed)

p_categorical<- plot_overlap(results_overlap_categorical)

p_smoking<- plot_overlap(results_overlap_smoking)



# plot directionality of ancestry-DMPs

library(dplyr)
ancestry_DML <- readRDS(paste0(project_path, "/Tissues/Ancestry_DML_signif.rds"))
ancestry_DML$tissue <- sub("\\..*", "", rownames(ancestry_DML))
ancestry_DML$direction <- ifelse(ancestry_DML$logFC > 0, "EA", "AA")
ancestry_DML$tissue <- factor(ancestry_DML$tissue, levels=rev(tissues))

ancestry_count <- ancestry_DML %>% group_by(tissue, direction) %>% count()
eur <- ancestry_DML %>% filter(logFC > 0) %>% group_by(tissue, direction) %>% count()
afr <- ancestry_DML %>% filter(logFC < 0) %>% group_by(tissue, direction) %>% count()


cols_ancestry <- c('EA'='#F0AE21','AA'='#F9DE8B')

direction_plot <- ggplot(ancestry_count, aes(x = tissue, y = n, fill = direction)) +
  geom_col(data = eur, aes(x = tissue, y = n, fill = direction), alpha = 1, position = "dodge") +
  geom_col(data = afr, aes(x = tissue, y = -n, fill = direction), alpha = 1, position = "dodge") +
  geom_hline(yintercept = 0, linetype = "dashed", linewidth = 0.3) +
  ylab("Number of CpGs") +xlab("") +coord_flip()+geom_text(data = afr, aes(label = n, x = tissue, y = -n), vjust = 0.5, hjust = 1, size = 3, position = position_dodge(width = 1)) +
  geom_text(data = eur, aes(label = n, x = tissue, y = n), hjust = -0, size = 3, position = position_dodge(width = 1)) +
  scale_fill_manual(values=cols_ancestry) +
  scale_y_continuous(labels = abs) + theme(legend.title = element_blank(),
                                            axis.text.x = element_text(colour="black", size=12),
                                            axis.text.y = element_text(colour="black", size=12),
                                            legend.text = element_text(colour="black", size=12),
                                            axis.title.x = element_text(size=12),
                                            legend.spacing.y = unit(-0.05, "cm"),
                                            panel.grid.major = element_blank(),
                                            panel.grid.minor = element_blank(),
                                            panel.border = element_rect(colour = "black", linewidth=1), legend.position="top")

pdf('~/cluster//Projects/GTEx_v8/Methylation/Plots/Ancestry_DMP_direction.pdf', 
    width = 4.68, height = 3.85)
print(direction_plot)
dev.off()


# plot smoking status per ancesrty ----
smoking_info <- read.table(paste0(project_path, "Donor_IDs_with_smoking_status.txt"), header = T)
colnames(smoking_info) <- c("SUBJID", "SmokerStatus", "Smoking")
tissues <- ctissues <- c("Lung", "ColonTransverse", "Ovary", "Prostate", "BreastMammaryTissue", "MuscleSkeletal", "KidneyCortex", "Testis", "WholeBlood")
mdata <- do.call(rbind.data.frame, lapply(tissues, function(tissue) readRDS(paste0(project_path, 'Tissues/', tissue, "/metadata.rds"))[,c("SUBJID","Ancestry")]))
mdata <- mdata[!duplicated(mdata$SUBJID),]
metadata <- merge(mdata, smoking_info, by='SUBJID')


ggplot(metadata[metadata$Ancestry!="AMR",], aes(x=Ancestry, alpha=SmokerStatus))+geom_bar(stat="count", position="fill")+scale_alpha_manual(values=c(0.4, 0.7, 1))+
  theme(legend.title = element_blank(),
        axis.text.x = element_text(colour="black", size=12),
        axis.text.y = element_text(colour="black", size=12),
        legend.text = element_text(colour="black", size=12),
        axis.title.x = element_text(size=12),
        legend.spacing.y = unit(-0.05, "cm"),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.border = element_rect(colour = "black", linewidth=1), legend.position="top", panel.background=element_blank())+theme(legend.position = "right")+ylab("Proportion")


# remove AMR and keep AFR/EUR
df <- metadata[metadata$Ancestry %in% c("AFR","EUR"), ]
# define smokers (adjust if your labels differ)
df[df$SmokerStatus == "ex-smoker",]$SmokerStatus <- "non-smoker"
df$Smoker <- df$SmokerStatus == "smoker"
tab <- table(df$Ancestry, df$Smoker)
tab

prop.test(
  x = tab[, "TRUE"],
  n = rowSums(tab)
)

# now do it per chromhmm

library(data.table)
library(dplyr)
library(valr)
library(ggplot2)
project_path <- "~/cluster/Projects/GTEx_v8/Methylation"
tissues <- c("Lung", "ColonTransverse", "Ovary", "Prostate", "BreastMammaryTissue", "MuscleSkeletal", "KidneyCortex", "Testis", "WholeBlood")

ancestry_DML <- readRDS(paste0(project_path,"/Tissues/Ancestry_DML_signif.rds"))

ancestry_DML$tissue <- sub("\\..*", "", rownames(ancestry_DML))
ancestry_DML$direction <- ifelse(ancestry_DML$logFC > 0, "EA", "AA")


#Load CpG coordinates (EPIC annotation)
annotation <- read.csv(paste0(project_path,"/Data/GPL21145_MethylationEPIC_15073387_v-1-0_processed.csv"))

ann_bed <- annotation[
  !is.na(annotation$MAPINFO) & !is.na(annotation$CHR),] %>%
  dplyr::select(chrom=CHR, start=MAPINFO, end=MAPINFO, name=IlmnID) %>%
  distinct()

ann_bed$chrom <- paste0("chr", ann_bed$chrom)
ann_bed$start <- ann_bed$start - 1

#chromHMM
chromhmm_dir <- paste0(project_path,"/ChromHMM/")
files <- list.files('/gpfs/projects/bsc83/Projects/GTEx_v8/Methylation/Data/EpiMap/', pattern='.bed.gz',full.names=T)

names_chrom <- c( "Lung","ColonTransverse","Ovary","Prostate",  "Breast","MuscleSkeletal","KidneyCortex","Testis","PBMC")

chromhmm <- lapply(names_chrom, function(tis){
  f <- files[grep(tis, files)]
  read.delim(f, header=FALSE, sep="\t")
  
})

names(chromhmm) <- names_chrom

# Intersect CpGs with ChromHMM regions

chromhmm_cpgs <- lapply(names_chrom, function(tis){
  chrom_df <- chromhmm[[tis]][,1:4]
  colnames(chrom_df) <- c("chrom","start","end","region")
  bed_intersect( ann_bed,chrom_df,suffix=c("_ann","_chromhmm")) %>% dplyr::select(cpg=name_ann, chromHMM=region_chromhmm) %>% distinct()
  
})

names(chromhmm_cpgs) <- names_chrom


# Assign ChromHMM state to ancestry DMP CpGs
ancestry_DML$cpg <- sub(".*\\.", "", rownames(ancestry_DML))
ancestry_DML$chromHMM <- NA

for(tis in names_chrom){
  idx <- ancestry_DML$tissue == tis
  ancestry_DML$chromHMM[idx] <- chromhmm_cpgs[[tis]]$chromHMM[match( ancestry_DML$cpg[idx],chromhmm_cpgs[[tis]]$cpg)  ]
}

# remove CpGs without ChromHMM annotation
ancestry_DML <- ancestry_DML[!is.na(ancestry_DML$chromHMM),]
saveRDS(ancestry_DML, paste0(project_path,"/Tissues/Ancestry_DML_signif_chromhmm.rds"))


# 6. Count CpGs per tissue / direction / ChromHMM
ancestry_DML_chromhmm <- readRDS(paste0(project_path,"/Tissues/Ancestry_DML_signif_chromhmm.rds"))
ancestry_DML_chromhmm$chromHMM[ancestry_DML_chromhmm$chromHMM %in% c("TssFlnkD", "TssFlnk", "TssFlnkU","TssA")] <- "TSS"
ancestry_DML_chromhmm$chromHMM[ancestry_DML_chromhmm$chromHMM %in% c("EnhA2", "EnhA1","EnhWk","EnhG1", "EnhG2")] <- "Enh"
ancestry_DML_chromhmm$chromHMM[ancestry_DML_chromhmm$chromHMM %in% c("ReprPCWk","ReprPC")] <- "ReprPC"
ancestry_DML_chromhmm$chromHMM[ancestry_DML_chromhmm$chromHMM %in% c("TxWk","Tx")] <- "Tx"
ancestry_count <- ancestry_DML_chromhmm %>% group_by(tissue, direction, chromHMM) %>%summarise(n = n(), .groups="drop")
ancestry_count$tissue <- factor(ancestry_count$tissue, levels = names_chrom)
eur <- ancestry_count %>% filter(direction=="EA")
afr <- ancestry_count %>% filter(direction=="AA")


cols_ancestry <- c("EA"="#F0AE21","AA"="#F9DE8B")

direction_plot <- ggplot() +geom_col( data=eur, aes(x=chromHMM, y=n, fill=direction), position="dodge") + 
  geom_col( data=afr, aes(x=chromHMM, y=-n, fill=direction), position="dodge") +
  geom_hline(yintercept=0, linetype="dashed", linewidth=0.3) +
  facet_wrap(~tissue, scales="free_y",nrow=3) + coord_flip() + scale_fill_manual(values=cols_ancestry) +  scale_y_continuous(labels=abs) +
  ylab("Number of CpGs") +  xlab("") +  theme_bw() +
  theme(legend.title=element_blank(),
    legend.position="top",
    panel.grid.major=element_blank(),
    panel.grid.minor=element_blank(),
    panel.border=element_rect(colour="black")
  )


pdf(paste0(project_path,"/Plots/Ancestry_DMP_direction_ChromHMM.pdf"),width=10,height=6)
print(direction_plot)
dev.off()


# plot smoking status
df <- metadata[metadata$Ancestry %in% c("AFR","EUR"), ]