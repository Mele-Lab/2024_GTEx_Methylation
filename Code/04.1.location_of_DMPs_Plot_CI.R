#!/usr/bin/env Rscript
# @Author: Winona Oliveros Diez
# @E-mail: winn95@gmail.com
# @Description: Plot fisher enrichments of DMPs per chromatin state and tissue
# @software version: R=4.2.2


library(tidyr)

####### Jose's code to plot both number and enrichment one plot per tissue and traitt

read_data <- function(variables, data, tissue, trait){ #Function to prepare data to plot and compute adjusted p value
  
  odds_ratio <- lapply(variables, function(type) data[[tissue]][[trait]][[type]][['f']]$estimate)
  adj.P.Val <- p.adjust(sapply(variables, function(type) data[[tissue]][[trait]][[type]][['f']]$p.value), method = "BH")
  CI_down <- lapply(variables, function(type) data[[tissue]][[trait]][[type]][['f']]$conf.int[1])
  CI_up <- lapply(variables, function(type) data[[tissue]][[trait]][[type]][['f']]$conf.int[2])
  sample_size <- lapply(variables, function(type) data[[tissue]][[trait]][[type]][['m']])
  
  
  names(odds_ratio) <- variables
  names(adj.P.Val) <- variables
  names(CI_down) <- variables
  names(CI_up) <- variables
  names(sample_size) <- variables
  
  odds_ratio_df <- as.data.frame(unlist(odds_ratio))
  odds_ratio_df$label <- variables
  odds_ratio_df$type <- deparse(substitute(data)) #Either hypo or hyper
  colnames(odds_ratio_df) <- c('oddsRatio', 'region','type')
  
  adj.P.Val_df <- as.data.frame(unlist(adj.P.Val))
  adj.P.Val_df$label <- variables
  adj.P.Val_df$type <- deparse(substitute(data))
  colnames(adj.P.Val_df) <- c('adjPvalue','region','type')
  
  CI_down_df <- as.data.frame(unlist(CI_down))
  CI_down_df$label <- variables
  CI_down_df$type <- deparse(substitute(data))
  colnames(CI_down_df) <- c('CI_down','region','type')
  
  CI_up_df <- as.data.frame(unlist(CI_up))
  CI_up_df$label <- variables
  CI_up_df$type <- deparse(substitute(data))
  colnames(CI_up_df) <- c('CI_up','region','type')
  
  sample_size_df <- as.data.frame(unlist(sample_size))
  sample_size_df$label <- variables
  sample_size_df$type <- deparse(substitute(data))
  colnames(sample_size_df) <- c('sample_size','region','type')
  
  all <- Reduce(function(x, y) merge(x, y, all=TRUE), list(odds_ratio_df, adj.P.Val_df, CI_down_df, CI_up_df, sample_size_df))
  head(all)
  all$sig <- 'not Sig'
  all$sig[all$adjPvalue<0.05] <- 'Sig'
  all <- all[,c("region","oddsRatio","adjPvalue","CI_down","CI_up","sig","type", "sample_size")]
  return(all)
}


hypo <- readRDS('~/marenostrum/Projects/GTEx_v8/Methylation/Tissues/enrichment_chromhmm_hypo_batch_CI.simple.continous.rds')
hyper <- readRDS('~/marenostrum/Projects/GTEx_v8/Methylation/Tissues/enrichment_chromhmm_hyper_batch_CI.simple.continous.rds')

ovary_hyper <- readRDS('~/marenostrum/Projects/GTEx_v8/Methylation/Tissues/Ovary.filt.enrichment_chromhmm_hyper.continous.rds')
ovary_hypo <- readRDS('~/marenostrum/Projects/GTEx_v8/Methylation/Tissues/Ovary.filt.enrichment_chromhmm_hypo.continous.rds')

colon_hyper <- readRDS('~/marenostrum/Projects/GTEx_v8/Methylation/Tissues/Colon.peer.enrichment_chromhmm_hyper.continous.rds')
colon_hypo <- readRDS('~/marenostrum/Projects/GTEx_v8/Methylation/Tissues/Colon.peers.enrichment_chromhmm_hypo.continous.rds')

hyper <- list()
hypo <- list()

hyper[['ColonTransverse']][['AGE']] <- colon_hyper
hypo[['ColonTransverse']][['AGE']] <- colon_hypo

colors_traits <- list('AGE'=c('#3D7CD0','#B4D6F6'),
                      'SEX2'=c('#3B734E','#89AA94'),
                      'EURv1'=c('#F0AE21','#F9DE8B'))

sex_tissues <- c('Ovary','Prostate','Testis')
#for (tissue in names(hypo)) {
for (tissue in c('ColonTransverse')) {
  #for (trait in names(hypo$Lung)) {
  for (trait in c('AGE')) {
    if (tissue %in% sex_tissues & trait == "SEX2") {
      print(NA)
    } else {
      hypo_d <- read_data(c("Enh","EnhBiv","Het","Quies","ReprPC","TSS","TssBiv","Tx","ZNF/Rpts"), hypo, tissue, trait)
      hyper_d <- read_data(c("Enh","EnhBiv","Het","Quies","ReprPC","TSS","TssBiv","Tx","ZNF/Rpts"), hyper, tissue, trait)
      hyper_hypo <- rbind(hypo_d, hyper_d)
      hyper_hypo$sig[hyper_hypo$sig =="Sig"] <- "FDR < 0.05"
      hyper_hypo$sig[hyper_hypo$sig =="not Sig"] <- "FDR >= 0.05"
      hyper_hypo$sig <- factor(hyper_hypo$sig, levels=c("FDR >= 0.05", "FDR < 0.05"))
      hyper_hypo$region <- factor(hyper_hypo$region, levels=rev(c("Enh","EnhBiv","Het","Quies","ReprPC","TSS","TssBiv","Tx","ZNF/Rpts")))
      hyper_hypo$type[hyper_hypo$type =="hypo"] <- "Hypomethylation"
      hyper_hypo$type[hyper_hypo$type =="hyper"] <- "Hypermethylation"
      g <- ggplot(hyper_hypo, aes(x=log2(oddsRatio), y=region, colour=type, alpha=sig)) +
        geom_errorbar(aes(xmin=log2(CI_down), xmax=log2(CI_up)), width=.3) +
        geom_vline(xintercept = 0) +
        #xlim(0,20) + #Only for Lung to show the 0
        geom_point(size=3) + ylab('') + theme_bw() +
        scale_colour_manual(values=colors_traits[[trait]]) +
        xlab("log2(Odds ratio)") +
        scale_alpha_discrete(range = c(0.4, 1), drop = FALSE) +
        theme(legend.title = element_blank(),
              axis.text.x = element_text(colour="black", size=13),
              axis.text.y = element_text(colour="black", size=14),
              legend.text = element_text(colour="black", size=13),
              axis.title.x = element_text(size=16),
              legend.spacing.y = unit(-0.05, "cm"),
              panel.grid.major = element_blank(),
              panel.grid.minor = element_blank(),
              panel.border = element_rect(colour = "black", linewidth=1)) +
        scale_y_discrete(breaks=c("Enh","EnhBiv","Het","Quies","ReprPC","TSS","TssBiv","Tx","ZNF/Rpts"),
                         labels=c("Enhancer","Enhancer Bivalent","Heterochromatin","Quiescent","Repressed Polycomb","TSS","TSS Bivalent","Transcription","ZNF & Repeats"))# + xlim(0, 3)

      #Plot sample sizes:
      
      g2 <- ggplot(hyper_hypo) + geom_col(aes(sample_size, region, fill=type), width = 0.6) +
        theme_classic() + xlab("Number of DMPs") + ylab("") +
        scale_fill_manual(values=colors_traits[[trait]]) +
        theme(legend.position = "none",
              axis.text.x = element_text(colour="black", size=13),
              axis.text.y=element_blank(),  #remove y axis labels,
              axis.title.x = element_text(size=16)) +
        scale_x_continuous(n.breaks=3)

      p <- ggarrange(g, g2, labels = c("A", "B"),
                     common.legend = TRUE, legend = "right", widths = c(0.8,0.3))
      pdf(file = paste0("~/marenostrum/Projects/GTEx_v8/Methylation/Plots/chromhmm/enrichment_", tissue,'_',trait,".v2.filt.pdf"), w = 8, h = 4)
      print(p)
      dev.off()
    }
  }
}


####### NOW plot per state and not per Tissue #######

read_data <- function(variables, data, type, trait){ #Function to prepare data to plot and compute adjusted p value
  
  odds_ratio <- lapply(variables, function(tissue) data[[tissue]][[trait]][[type]][['f']]$estimate)
  adj.P.Val <- p.adjust(sapply(variables, function(tissue) data[[tissue]][[trait]][[type]][['f']]$p.value), method = "BH")
  CI_down <- lapply(variables, function(tissue) data[[tissue]][[trait]][[type]][['f']]$conf.int[1])
  CI_up <- lapply(variables, function(tissue) data[[tissue]][[trait]][[type]][['f']]$conf.int[2])
  sample_size <- lapply(variables, function(tissue) data[[tissue]][[trait]][[type]][['m']])
  
  
  names(odds_ratio) <- variables
  names(adj.P.Val) <- variables
  names(CI_down) <- variables
  names(CI_up) <- variables
  names(sample_size) <- variables
  
  odds_ratio_df <- as.data.frame(unlist(odds_ratio))
  odds_ratio_df$label <- variables
  odds_ratio_df$type <- deparse(substitute(data)) #Either hypo or hyper
  colnames(odds_ratio_df) <- c('oddsRatio', 'tissue','type')
  
  adj.P.Val_df <- as.data.frame(unlist(adj.P.Val))
  adj.P.Val_df$label <- variables
  adj.P.Val_df$type <- deparse(substitute(data))
  colnames(adj.P.Val_df) <- c('adjPvalue','tissue','type')
  
  CI_down_df <- as.data.frame(unlist(CI_down))
  CI_down_df$label <- variables
  CI_down_df$type <- deparse(substitute(data))
  colnames(CI_down_df) <- c('CI_down','tissue','type')
  
  CI_up_df <- as.data.frame(unlist(CI_up))
  CI_up_df$label <- variables
  CI_up_df$type <- deparse(substitute(data))
  colnames(CI_up_df) <- c('CI_up','tissue','type')
  
  sample_size_df <- as.data.frame(unlist(sample_size))
  sample_size_df$label <- variables
  sample_size_df$type <- deparse(substitute(data))
  colnames(sample_size_df) <- c('sample_size','tissue','type')
  
  all <- Reduce(function(x, y) merge(x, y, all=TRUE), list(odds_ratio_df, adj.P.Val_df, CI_down_df, CI_up_df, sample_size_df))
  head(all)
  all$sig <- 'not Sig'
  all$sig[all$adjPvalue<0.05] <- 'Sig'
  all <- all[,c("tissue","oddsRatio","adjPvalue","CI_down","CI_up","sig","type", "sample_size")]
  return(all)
}


hypo <- readRDS('~/marenostrum/Projects/GTEx_v8/Methylation/Tissues/enrichment_chromhmm_hypo_batch_CI.simple.continous.rds')
hyper <- readRDS('~/marenostrum/Projects/GTEx_v8/Methylation/Tissues/enrichment_chromhmm_hyper_batch_CI.simple.continous.rds')

colors_traits <- list('AGE'=c('#3D7CD0','#B4D6F6'),
                      'SEX2'=c('#3B734E','#89AA94'),
                      'EURv1'=c('#F0AE21','#F9DE8B'))

tissue_info <- readRDS(paste0('~/marenostrum/', "Projects/GTEx_v8/Jose/00_Data/Tissue_info_whole.rds"))

tissues <- c("Lung", "ColonTransverse", "Ovary", "Prostate", "BreastMammaryTissue", "KidneyCortex", "Testis", "WholeBlood","MuscleSkeletal")

library(ggpubr)
sex_tissues <- c('Ovary','Prostate','Testis')
for (type in c("Enh","EnhBiv","Het","Quies","ReprPC","TSS","TssBiv","Tx","ZNF/Rpts")) {
  for (trait in c('AGE','SEX2','EURv1')) {
    if (trait == 'SEX2') {
      tissues_to_test <- tissues[!tissues %in% sex_tissues]
    } 
    if (trait == 'AGE') {
      tissues_to_test <- tissues[tissues != 'MuscleSkeletal']
    }
    if (trait == 'EURv1') {
      tissues_to_test <- tissues
    }
      hypo_d <- read_data(tissues_to_test, hypo, type, trait)
      hyper_d <- read_data(tissues_to_test, hyper, type, trait)
      hyper_hypo <- rbind(hypo_d, hyper_d)
      hyper_hypo$sig[hyper_hypo$sig =="Sig"] <- "FDR < 0.05"
      hyper_hypo$sig[hyper_hypo$sig =="not Sig"] <- "FDR >= 0.05"
      hyper_hypo$sig <- factor(hyper_hypo$sig, levels=c("FDR >= 0.05", "FDR < 0.05"))
      hyper_hypo$tissue <- factor(hyper_hypo$tissue, levels=tissues_to_test)
      hyper_hypo$type[hyper_hypo$type =="hypo"] <- "Hypomethylation"
      hyper_hypo$type[hyper_hypo$type =="hyper"] <- "Hypermethylation"
      g <- ggplot(hyper_hypo, aes(x=log2(oddsRatio), y=tissue, colour=type, alpha=sig)) +
        geom_errorbar(aes(xmin=log2(CI_down), xmax=log2(CI_up)), width=.3) +
        geom_vline(xintercept = 0) +
        #xlim(0,20) + #Only for Lung to show the 0
        geom_point(size=3) + ylab('') + theme_bw() +
        scale_colour_manual(values=colors_traits[[trait]]) +
        xlab("log2(Odds ratio)") +
        scale_alpha_discrete(range = c(0.4, 1), drop = FALSE) +
        theme(legend.title = element_blank(),
              axis.text.x = element_text(colour="black", size=13),
              axis.text.y = element_text(colour="black", size=14),
              legend.text = element_text(colour="black", size=13),
              axis.title.x = element_text(size=16),
              legend.spacing.y = unit(-0.05, "cm"),
              panel.grid.major = element_blank(),
              panel.grid.minor = element_blank(),
              panel.border = element_rect(colour = "black", linewidth=1)) #+

      #Plot sample sizes:
      
      g2 <- ggplot(hyper_hypo) + geom_col(aes(sample_size, tissue, fill=type), width = 0.6) +
        theme_classic() + xlab("Number of DMPs") + ylab("") +
        scale_fill_manual(values=colors_traits[[trait]]) +
        theme(legend.position = "none",
              axis.text.x = element_text(colour="black", size=13),
              axis.text.y=element_blank(),  #remove y axis labels,
              axis.title.x = element_text(size=16)) +
              scale_x_continuous(n.breaks=3)

      p <- ggarrange(g, g2, labels = c("A", "B"),
                     common.legend = TRUE, legend = "right", widths = c(0.8,0.3))
      pdf(file = paste0("~/marenostrum/Projects/GTEx_v8/Methylation/Plots/chromhmm/enrichment", gsub('\\/','_',type),'_',trait,".v2.pdf"), w = 8, h = 4)
      print(p)
      dev.off()
    }
}


###### enrichment shared positions #####
read_data <- function(variables, data, trait){ #Function to prepare data to plot and compute adjusted p value
  
  odds_ratio <- lapply(variables, function(type) data[[trait]][[type]][['f']]$estimate)
  adj.P.Val <- p.adjust(sapply(variables, function(type) data[[trait]][[type]][['f']]$p.value), method = "BH")
  CI_down <- lapply(variables, function(type) data[[trait]][[type]][['f']]$conf.int[1])
  CI_up <- lapply(variables, function(type) data[[trait]][[type]][['f']]$conf.int[2])
  sample_size <- lapply(variables, function(type) data[[trait]][[type]][['m']])
  
  
  names(odds_ratio) <- variables
  names(adj.P.Val) <- variables
  names(CI_down) <- variables
  names(CI_up) <- variables
  names(sample_size) <- variables
  
  odds_ratio_df <- as.data.frame(unlist(odds_ratio))
  odds_ratio_df$label <- variables
  odds_ratio_df$type <- deparse(substitute(data)) #Either hypo or hyper
  colnames(odds_ratio_df) <- c('oddsRatio', 'region','type')
  
  adj.P.Val_df <- as.data.frame(unlist(adj.P.Val))
  adj.P.Val_df$label <- variables
  adj.P.Val_df$type <- deparse(substitute(data))
  colnames(adj.P.Val_df) <- c('adjPvalue','region','type')
  
  CI_down_df <- as.data.frame(unlist(CI_down))
  CI_down_df$label <- variables
  CI_down_df$type <- deparse(substitute(data))
  colnames(CI_down_df) <- c('CI_down','region','type')
  
  CI_up_df <- as.data.frame(unlist(CI_up))
  CI_up_df$label <- variables
  CI_up_df$type <- deparse(substitute(data))
  colnames(CI_up_df) <- c('CI_up','region','type')
  
  sample_size_df <- as.data.frame(unlist(sample_size))
  sample_size_df$label <- variables
  sample_size_df$type <- deparse(substitute(data))
  colnames(sample_size_df) <- c('sample_size','region','type')
  
  all <- Reduce(function(x, y) merge(x, y, all=TRUE), list(odds_ratio_df, adj.P.Val_df, CI_down_df, CI_up_df, sample_size_df))
  head(all)
  all$sig <- 'not Sig'
  all$sig[all$adjPvalue<0.05] <- 'Sig'
  all <- all[,c("region","oddsRatio","adjPvalue","CI_down","CI_up","sig","type", "sample_size")]
  return(all)
}

##### plot shared positions
# hypo <- readRDS('~/marenostrum/Projects/GTEx_v8/Methylation/Tissues/enrichment_chromhmm_hypo_shared_CI.continous.2.rds')
# hyper <- readRDS('~/marenostrum/Projects/GTEx_v8/Methylation/Tissues/enrichment_chromhmm_hyper_shared_CI.continous.2.rds')

hypo <- readRDS(paste0(basepath,'/Projects/GTEx_v8/Methylation/Tissues/enrichment_chromhmm_hypo_shared_CI.continous.2.rds'))
hyper <- readRDS(paste0(basepath,'/Projects/GTEx_v8/Methylation/Tissues/enrichment_chromhmm_hyper_shared_CI.continous.2.rds'))


colors_traits <- list('Age'=c('#3D7CD0','#B4D6F6'),
                      'Sex'=c('#3B734E','#89AA94'),
                      'Ancestry'=c('#F0AE21','#F9DE8B'))

library(ggpubr)
  for (trait in c('Ancestry','Sex','Age')) {

    hypo_d <- read_data(c("Enh","EnhBiv","Het","Quies","ReprPC","TSS","TssBiv","Tx","ZNF/Rpts"), hypo, trait)
    hyper_d <- read_data(c("Enh","EnhBiv","Het","Quies","ReprPC","TSS","TssBiv","Tx","ZNF/Rpts"), hyper, trait)
    hyper_hypo <- rbind(hypo_d, hyper_d)
    hyper_hypo$sig[hyper_hypo$sig =="Sig"] <- "FDR < 0.05"
    hyper_hypo$sig[hyper_hypo$sig =="not Sig"] <- "FDR >= 0.05"
    hyper_hypo$sig <- factor(hyper_hypo$sig, levels=c("FDR >= 0.05", "FDR < 0.05"))
    hyper_hypo$region <- factor(hyper_hypo$region, levels=rev(c("Enh","EnhBiv","Het","Quies","ReprPC","TSS","TssBiv","Tx","ZNF/Rpts")))
    hyper_hypo$type[hyper_hypo$type =="hypo"] <- "Hypomethylation"
    hyper_hypo$type[hyper_hypo$type =="hyper"] <- "Hypermethylation"
    g <- ggplot(hyper_hypo, aes(x=log2(oddsRatio), y=region, colour=type, alpha=sig)) +
      geom_errorbar(aes(xmin=log2(CI_down), xmax=log2(CI_up)), width=.3) +
      geom_vline(xintercept = 0) +
      #xlim(0,20) + #Only for Lung to show the 0
      geom_point(size=3) + ylab('') + theme_bw() +
      scale_colour_manual(values=colors_traits[[trait]]) +
      xlab("log2(Odds ratio)") +
      scale_alpha_discrete(range = c(0.4, 1), drop = FALSE) +
      theme(legend.title = element_blank(),
            axis.text.x = element_text(colour="black", size=13),
            axis.text.y = element_text(colour="black", size=14),
            legend.text = element_text(colour="black", size=13),
            axis.title.x = element_text(size=16),
            legend.spacing.y = unit(-0.05, "cm"),
            panel.grid.major = element_blank(),
            panel.grid.minor = element_blank(),
            panel.border = element_rect(colour = "black", linewidth=1)) +
      scale_y_discrete(breaks=c("Enh","EnhBiv","Het","Quies","ReprPC","TSS","TssBiv","Tx","ZNF/Rpts"),
                       labels=c("Enhancer","Enhancer Bivalent","Heterochromatin","Quiescent","Repressed Polycomb","TSS","TSS Bivalent","Transcription","ZNF & Repeats"))# + xlim(0, 3)

    #Plot sample sizes:
    
    g2 <- ggplot(hyper_hypo) + geom_col(aes(sample_size, region, fill=type), width = 0.6) +
      theme_classic() + xlab("Number of DMPs") + ylab("") +
      scale_fill_manual(values=colors_traits[[trait]]) +
      theme(legend.position = "none",
            axis.text.x = element_text(colour="black", size=13),
            axis.text.y=element_blank(),
            axis.title.x = element_text(size=16)) +
      scale_y_discrete(breaks=c("Enh","EnhBiv","Het","Quies","ReprPC","TSS","TssBiv","Tx","ZNF/Rpts"),
                       labels=c("Enhancer","Enhancer Bivalent","Heterochromatin","Quiescent","Repressed Polycomb","TSS","TSS Bivalent","Transcription","ZNF & Repeats")) #+

    p <- ggarrange(g, g2, labels = c("A", "B"),
                   common.legend = TRUE, legend = "right", widths = c(0.8,0.3))
    pdf(file = paste0("~/marenostrum/Projects/GTEx_v8/Methylation/Plots/chromhmm/enrichment",'_',trait,"_shared.2.pdf"), w = 8, h = 4)
    print(p)
    dev.off()
  }


# plot enrichment for tissue-specific positions 
read_data_simple <- function(regions, fisher_list) {
  odds_ratio <- sapply(regions, function(r) fisher_list[[r]][["f"]]$estimate)
  pvals <- sapply(regions, function(r) fisher_list[[r]][["f"]]$p.value)
  adj.P.Val <- p.adjust(pvals, method="BH")
  CI_down <- sapply(regions, function(r) fisher_list[[r]][["f"]]$conf.int[1])
  CI_up <- sapply(regions, function(r) fisher_list[[r]][["f"]]$conf.int[2])
  sample_size <- sapply(regions, function(r) fisher_list[[r]][["m"]])
  
  dt <- data.frame(region=regions,
                   oddsRatio=odds_ratio,
                   adjPvalue=adj.P.Val,
                   CI_down=CI_down,
                   CI_up=CI_up,
                   sample_size=sample_size)
  dt$sig <- ifelse(dt$adjPvalue < 0.05, "FDR < 0.05", "FDR >= 0.05")
  dt$sig <- factor(dt$sig, levels=c("FDR >= 0.05", "FDR < 0.05"))
  dt$region <- factor(dt$region, levels=rev(regions))
  dt
}

library(ggplot2)
library(ggpubr)

families <- c('Enh','EnhBiv','Het','Quies','ReprPC','TSS','TssBiv','Tx','ZNF/Rpts')
dt_AA <- read_data_simple(families, fisher_results_nonshared[["AA"]])
dt_AA$type <- "AA"
dt_EA <- read_data_simple(families, fisher_results_nonshared[["EA"]])
dt_EA$type <- "EA"
dt <- rbind(dt_EA, dt_AA)

g <- ggplot(dt, aes(x=log2(oddsRatio), y=region, alpha=sig,color=type)) +
  geom_errorbar(aes(xmin=log2(CI_down), xmax=log2(CI_up)), width=.3) +
  geom_vline(xintercept = 0) +
  geom_point(size=3) + ylab('') + theme_bw() +
  scale_colour_manual(values=rev(colors_traits[["Ancestry"]])) +
  xlab("log2(Odds ratio)") +
  scale_alpha_discrete(range = c(0.4, 1), drop = FALSE) +
  theme(legend.title = element_blank(),
        axis.text.x = element_text(colour="black", size=13),
        axis.text.y = element_text(colour="black", size=14),
        legend.text = element_text(colour="black", size=13),
        axis.title.x = element_text(size=16),
        legend.spacing.y = unit(-0.05, "cm"),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.border = element_rect(colour = "black", linewidth=1)) +
  scale_y_discrete(breaks=c("Enh","EnhBiv","Het","Quies","ReprPC","TSS","TssBiv","Tx","ZNF/Rpts"),
                   labels=c("Enhancer","Enhancer Bivalent","Heterochromatin","Quiescent","Repressed Polycomb","TSS","TSS Bivalent","Transcription","ZNF & Repeats"))+
  ggtitle("Tissue-specific DMPs")# + xlim(0, 3)


g2 <- ggplot(dt) + geom_col(aes(sample_size, region, fill=type), width = 0.6) +
  theme_classic() + xlab("Number of DMPs") + ylab("") +
  scale_fill_manual(values=rev(colors_traits[["Ancestry"]])) +
  theme(legend.position = "none",
        axis.text.x = element_text(colour="black", size=13),
        axis.text.y=element_blank(),
        axis.title.x = element_text(size=16)) +
  scale_y_discrete(breaks=c("Enh","EnhBiv","Het","Quies","ReprPC","TSS","TssBiv","Tx","ZNF/Rpts"),
                   labels=c("Enhancer","Enhancer Bivalent","Heterochromatin","Quiescent","Repressed Polycomb","TSS","TSS Bivalent","Transcription","ZNF & Repeats")) #+

p <- ggarrange(g, g2, widths=c(0.8, 0.3))
pdf(paste0(basepath, "/Projects/GTEx_v8/Methylation/Plots/chromhmm/enrichment_nonshared.pdf"),
    w=8, h=4)
print(p)
dev.off()
