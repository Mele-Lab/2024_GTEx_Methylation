#!/usr/bin/env Rscript
# @Author: Jose Miguel Ramirez; Adapted by Winona Oliveros
# @E-mail: jose.ramirez1@bsc.es
# @Description: Code to parse correlations for later plotting
# @software version: R=4.2.2

#Set path 
# setwd(system("pwd", intern = T)) #If in linux
# setwd(dirname(rstudioapi::getActiveDocumentContext()$path)) #If using Rstudio
# setwd("..")

# Tissues ----
### read methylation results ####
first_dir <- "/gpfs/projects/bsc83/"
project_path <- paste0(first_dir, "Projects/GTEx_v8/Methylation/")

tissues <- c("Lung", "ColonTransverse", "Ovary", "Prostate", "BreastMammaryTissue", "MuscleSkeletal", "KidneyCortex", "Testis", "WholeBlood")
names <- c("Age", "Ancestry", "BMI", "Sex")

tissue_info <- readRDS(paste0(first_dir, "Projects/GTEx_v8/Methylation/Data/Tissue_info_whole.rds"))

tissue_info <- tissue_info[!grepl("BreastMammaryTissue_", tissue_info$tissue_ID),]
tissue_info <- tissue_info[tissue_info$tissue_ID %in% tissues,]

sex_tissues <- c('Ovary','Prostate','Testis')

tissues <- tissue_info$tissue_ID

### read Correlations
get_corr <- function(tissue, trait){
  if(tissue %in% sex_tissues & trait == "SEX2"){
    NA
  }else{
    model <- readRDS(paste0(project_path, "Tissues/",tissue, "/",trait,'_Correlations_probes_genes_DEG_DMP.padjusted_ancestry_c.rds'))
    #rownames(model[[trait]][model[[trait]]$adj.P.Val<0.05,])
    model
  }
}
DMPs_Res <- lapply(c('EURv1','SEX2','AGE','BMI'), function(trait) lapply(tissues, function(tissue) get_corr(tissue, trait)))
names(DMPs_Res) <- c("Ancestry", "Sex", "Age", "BMI")
for(trait in c("Ancestry", "Sex", "Age", "BMI")){names(DMPs_Res[[trait]]) <- tissues}

get_sig <- function(tissue, trait){
  print(paste0(trait,':',tissue))
  if(tissue %in% sex_tissues & trait == "Sex"){
    NA
  }else{
    DMPs_Res[[trait]][[tissue]] <- DMPs_Res[[trait]][[tissue]][!is.na(DMPs_Res[[trait]][[tissue]]$gene),]
    table(DMPs_Res[[trait]][[tissue]]$p.adj<0.05)
  }
}

lapply(c("Ancestry", "Sex", "Age", "BMI"), function(trait) lapply(tissues, function(tissue) get_sig(tissue, trait)))

get_signif <- function(tissue, trait){
  if(tissue %in% sex_tissues & trait == "Sex"){
    NA
  }else{
    DMPs_Res[[trait]][[tissue]][DMPs_Res[[trait]][[tissue]]$p.adj<0.05,]
  }
}
signif <- lapply(c("Ancestry", "Sex", "Age", "BMI"), function(trait) lapply(tissues, function(tissue) get_signif(tissue, trait)))
names(signif) <- c("Ancestry", "Sex", "Age", "BMI")
for(trait in c("Ancestry", "Sex", "Age", "BMI")){names(signif[[trait]]) <- tissues}

#Plot rho values:
plot <- function(tissue, trait){
  if(tissue %in% sex_tissues & trait == "Sex"){
    NA
  }else{
    ggplot(signif[[trait]][[tissue]], aes(class, abs(cor))) + geom_violin() + 
      geom_boxplot(outlier.shape = NA,
                   notch = T,
                   width = 0.25) + theme_bw() + xlab("") + ylab("abs(rho)") +
  ggtitle(paste0(tissue,':',trait))
  }
}

library(ggplot2)
lapply(c("Ancestry", "Sex", "Age", "BMI"), function(trait) lapply(tissues, function(tissue) plot(tissue, trait)))

### read DEG ####

to_plot_2 <- c(1:6)
to_plot_1 <- c(1:6)
for (tissue in tissues) {
  for (trait in c("Ancestry", "Sex", "Age", "BMI")) {
    if(tissue %in% sex_tissues & trait == "Sex"){
      next
    } else {
      
        #Promoters
        pro_all <- DMPs_Res[[trait]][[tissue]][!is.na(DMPs_Res[[trait]][[tissue]]$gene) & DMPs_Res[[trait]][[tissue]]$class=="promoter",]
        #pro_deg <- pro_all[pro_all$gene %in% degs,] #DEG-probe pairs
        bg <- length(unique(pro_all$gene)) # DEGs that have at least one probe in promoters
        corr <- pro_all[pro_all$p.adj<0.05,]
        
        ### gene numbers 
        to_plot_2 <- rbind(to_plot_2,c(length(unique(corr[corr$cor>0,"gene"]))/bg, "Promoter", "Positive", tissue, trait,length(unique(corr[corr$cor>0,"gene"]))),
                           c(length(unique(corr[corr$cor<0,"gene"]))/bg, "Promoter", "Negative", tissue, trait,length(unique(corr[corr$cor<0,"gene"]))))

        bg <- length(unique(pro_all$probe)) # DMPs that are associated to a gene
        ### probe numbers 
        to_plot_1 <- rbind(to_plot_1,c(length(unique(corr[corr$cor>0,"probe"]))/bg, "Promoter", "Positive", tissue, trait,length(unique(corr[corr$cor>0,"probe"]))),
                           c(length(unique(corr[corr$cor<0,"probe"]))/bg, "Promoter", "Negative", tissue, trait,length(unique(corr[corr$cor<0,"probe"]))))
        
        #Enhancers
        pro_all <- DMPs_Res[[trait]][[tissue]][!is.na(DMPs_Res[[trait]][[tissue]]$gene) & DMPs_Res[[trait]][[tissue]]$class=="enhancer",]
        #pro_deg <- pro_all[pro_all$gene %in% degs,] #DEG-probe pairs
        bg <- length(unique(pro_all$gene)) # DEGs that have at least one probe in promoters
        corr <- pro_all[pro_all$p.adj<0.05,]
        
        ### gene numbers 
        to_plot_2 <- rbind(to_plot_2,c(length(unique(corr[corr$cor>0,"gene"]))/bg, "Enhancer", "Positive", tissue, trait,length(unique(corr[corr$cor>0,"gene"]))),
                           c(length(unique(corr[corr$cor<0,"gene"]))/bg, "Enhancer", "Negative", tissue, trait,length(unique(corr[corr$cor<0,"gene"]))))
        
        bg <- length(unique(pro_all$probe)) # DMPs that are associated to a gene
        ### probe numbers 
        to_plot_1 <- rbind(to_plot_1,c(length(unique(corr[corr$cor>0,"probe"]))/bg, "Enhancer", "Positive", tissue, trait,length(unique(corr[corr$cor>0,"probe"]))),
                           c(length(unique(corr[corr$cor<0,"probe"]))/bg, "Enhancer", "Negative", tissue, trait,length(unique(corr[corr$cor<0,"probe"]))))
        
        #Gene Body
        pro_all <- DMPs_Res[[trait]][[tissue]][!is.na(DMPs_Res[[trait]][[tissue]]$gene) & DMPs_Res[[trait]][[tissue]]$class=="gene_body",]
        #pro_deg <- pro_all[pro_all$gene %in% degs,] #DEG-probe pairs
        bg <- length(unique(pro_all$gene)) # DEGs that have at least one probe in promoters
        corr <- pro_all[pro_all$p.adj<0.05,]
        
        ### gene numbers 
        to_plot_2 <- rbind(to_plot_2,c(length(unique(corr[corr$cor>0,"gene"]))/bg, "Gene Body", "Positive", tissue, trait,length(unique(corr[corr$cor>0,"gene"]))),
                           c(length(unique(corr[corr$cor<0,"gene"]))/bg, "Gene Body", "Negative", tissue, trait,length(unique(corr[corr$cor<0,"gene"]))))
        
        bg <- length(unique(pro_all$probe)) # DMPs that are associated to a gene
        ### probe numbers 
        to_plot_1 <- rbind(to_plot_1,c(length(unique(corr[corr$cor>0,"probe"]))/bg, "Gene Body", "Positive", tissue, trait,length(unique(corr[corr$cor>0,"probe"]))),
                           c(length(unique(corr[corr$cor<0,"probe"]))/bg, "Gene Body", "Negative", tissue, trait,length(unique(corr[corr$cor<0,"probe"]))))
        
      }
    }
  }


saveRDS(to_plot_2, paste0(project_path, '/Data/correlations_to_plot_2_ancestry_continous.new.rds'))
saveRDS(to_plot_1, paste0(project_path, "/Data/correlations_to_plot_1_ancestry_continous.new.rds"))


