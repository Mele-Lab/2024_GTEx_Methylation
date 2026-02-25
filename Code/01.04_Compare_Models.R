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

ancestry_DML <- readRDS(paste0(project_path, "/Tissues/Ancestry_DML_signif.rds"))
ancestry_DML$tissue <- sub("\\..*", "", rownames(ancestry_DML))

ggplot()





