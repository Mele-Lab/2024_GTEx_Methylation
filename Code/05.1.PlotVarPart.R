### Plot varPart #####
library(variancePartition)
first_dir <- "~/"
#setwd(paste0(first_dir, "marenostrum/Projects/GTEx_v8/Methylation/"))
basepath <- "/Users/mariasopenar/cluster/"
setwd(paste0(basepath, "/Projects/GTEx_v8/Methylation/"))

# -------------- #
print(Sys.time())
#-------------- #

tissues <- list.dirs("Tissues/", full.names = F)[-1]
tissues <- tissues[-grep('Old',tissues)]
chuncks <- c(1:16)

#### reading varPart values ####
beta <- lapply(chuncks, function(chnk) readRDS(paste0('varPart/', chnk, "_chunck_var_part.rds")))
beta_df <- do.call("rbind",beta)

#### plot 
head(beta_df)
pdf(file = paste0('Plots/',"/VarPart_general.pdf"), w = 6, h = 4)
plotVarPart( sortCols(beta_df))
dev.off()

#### read DM results ####
results_DML <- lapply(tissues, function(tis) 
  readRDS(paste0("/Projects/GTEx_v8/Methylation/Tissues/",tis,"/DML_results_5_PEERs_continous.rds")))
names(results_DML) <- tissues

results_DML_all <- lapply(tissues, function(tis) 
  do.call("rbind",results_DML[[tis]]))
names(results_DML_all) <- tissues

## sig results 
DMP <- lapply(tissues, function(tis) 
  rownames(results_DML_all[[tis]][results_DML_all[[tis]]$adj.P.Val<0.05,]))
names(DMP) <- tissues

### all DMP - union ####
all_dmp <- unique(gsub('.*\\.','', (unlist(DMP))))

##plot
head(beta_df)
pdf(file = paste0('Plots/',"/VarPart_general.DMP.pdf"), w = 6, h = 4)
plotVarPart( sortCols(beta_df[all_dmp,]))
dev.off()

### per trait #### 
for (trait in c('EURv1','SEX2','AGE','BMI')) {
  cpgs <- unlist(DMP)[grep(trait, unlist(DMP))]
  names <- unique(gsub('.*\\.','', cpgs))
  
  pdf(file = paste0('Plots/',trait,"_VarPart_general.DMP.pdf"), w = 6, h = 4)
  print(plotVarPart( sortCols(beta_df[names,])))
  dev.off()
}

### merge info with prom/enh 
Sys.time()
#data_path <- "~/marenostrum/Projects/GTEx_v8/Methylation/Data/"
#annotation <- read.delim(paste0(data_path, "Methylation_Epic_gene_promoter_enhancer_processed.txt"), sep = '\t', header = T)
annotation <- read.delim(paste0(basepath, "Projects/GTEx_v8/Methylation/Data/Methylation_Epic_gene_promoter_enhancer_processed.txt"), sep = '\t', header = T)
Sys.time()

#Running a different analysis per promoter, enhancer and gene body
beta_df$Type <- 'Other'
beta_df$Type[rownames(beta_df) %in% annotation$IlmnID[annotation$Type=="Promoter_Associated"]] <- "Promoter_Associated"
beta_df$Type[rownames(beta_df) %in% annotation$IlmnID[annotation$Type=="Enhancer_Associated"]] <- "Enhancer_Associated"
beta_df$Type[rownames(beta_df) %in% annotation$IlmnID[annotation$Type=="Gene_Associated"]] <-"Gene_Associated"
table(beta_df$Type)

ggplot(beta_df, aes(x=Tissue, y=SUBJID, color=Type)) + 
  geom_point(alpha = 0.3) + 
  #scale_color_manual(values=c( "#95B46A","#F98948"))+
  #geom_point(data=expr_variation[expr_variation$Ribosomal == 'yes',], aes(x=Tissue, y=Individual), color='black') +
  #geom_text(aes(label=ifelse(Ribosomal == 'yes',as.character(gene_name),'')),hjust=0,vjust=0, color = 'black') +
  # geom_text_repel(
  #   data = subset(expr_variation, Ribosomal == 'yes'),
  #   aes(label = gene_name),
  #   size = 3,
  #   segment.color = 'grey50',
  #   color='black'
  # ) +
  theme_bw()

ggplot(beta_df, aes(x=Tissue, y=SUBJID) ) +
  geom_hex(bins = 70) +
  scale_fill_continuous(type = "viridis") + facet_grid(~ Type)+
  theme_bw()

library(dplyr)
library(purrr)
library(ggplot2)
library(cowplot)

pdf(file = paste0('Plots/',"/VarPart_per_loc.pdf"), w = 15, h = 3.5)
beta_df %>% 
  group_split(Type) %>% 
  map(
    ~ggplot(as.data.frame(.), aes(x=Tissue, y=SUBJID)) + 
      geom_hex(bins = 70) +
      scale_fill_continuous(type = "viridis") +
      theme_bw() +
      facet_grid(~ Type)
  ) %>% 
  plot_grid(plotlist = ., ncol = 4)
dev.off()

for (type in c('Promoter_Associated',"Enhancer_Associated",'Gene_Associated')) {
  plt <- ggplot(beta_df[beta_df$Type == type,], aes(x=Tissue, y=SUBJID, color=Type)) +
    geom_point() +
    #scale_color_manual(values=c( "#95B46A","#F98948"))+
    #geom_point(data=expr_variation[expr_variation$Ribosomal == 'yes',], aes(x=Tissue, y=Individual), color='black') +
    #geom_text(aes(label=ifelse(Ribosomal == 'yes',as.character(gene_name),'')),hjust=0,vjust=0, color = 'black') +
    # geom_text_repel(
    #   data = subset(expr_variation, Ribosomal == 'yes'),
    #   aes(label = gene_name),
    #   size = 3,
    #   segment.color = 'grey50',
    #   color='black'
    # ) +
    theme_bw()
  print(plt)
  # plt2 <-  ggplot(beta_df[beta_df$Type == type,], aes(x=Tissue, y=SUBJID) ) +
  #   stat_density_2d(aes(fill = ..level..), geom = "polygon") + 
  #   ggtitle(type)
  # print(plt2)
  
  plt2 <- ggplot(beta_df[beta_df$Type == type,], aes(x=Tissue, y=SUBJID) ) +
    geom_hex(bins = 70) +
    scale_fill_continuous(type = "viridis") +
    theme_bw()
  
  print(plt2)
}


#### merge data 
beta_df$IlmnID <- rownames(beta_df)
beta_df_info <- merge(beta_df, annotation, by='IlmnID', all.x=TRUE)

### enrichment genomic locations
shared_cpgs <- readRDS('~/marenostrum/Projects/GTEx_v8/Methylation/Data/chromHMM_shared_9tissues.rds')
beta_df$cpg <- rownames(beta_df)
betas_locations <- merge(beta_df, shared_cpgs, by.x = 'cpg', by.y = 'name_ann')
head(betas_locations)

library(reshape2)
betas_locations_m <- melt(data = betas_locations[,c("cpg", "SUBJID", "Tissue", "region_chromhmm")], id.vars = c('cpg','region_chromhmm'),
                          variable.name = 'Variable', value.name = 'Proportion')
head(betas_locations_m)

library(ggplot2)
# library(extrafont)
# font_import()
pdf(file = paste0('Plots/',"/VarPart_chromhmm.pdf"), w = 5, h = 3)
ggplot(betas_locations_m, aes(x = Proportion, y = region_chromhmm)) + 
  geom_violin(fill='grey')+ facet_grid(~ Variable)+
  theme_bw()
dev.off()


#### main plot overview varPart ####
## general plot
beta_df$cpg <- rownames(beta_df)

library(reshape2)
betas_locations_m <- melt(data = beta_df[,c("cpg", "SUBJID", "Tissue")], id.vars = c('cpg'),
                          variable.name = 'Variable', value.name = 'Proportion')
head(betas_locations_m)

library("dplyr")
mu <- betas_locations_m %>% 
  group_by(Variable) %>%
  summarise(grp.mean = mean(Proportion))
mu

pdf(file = paste0('Plots/',"/VarPart_all_main2_meth.pdf"), w = 4, h = 3)
ggplot(betas_locations_m, aes(x = Proportion, fill=Variable)) + 
  geom_density() +
  scale_fill_brewer(palette="Dark2") +
  geom_vline(data=mu, aes(xintercept=grp.mean, color=Variable),
             linetype="dashed") +
  geom_text(
    size    = 5,
    data    = mu,
    mapping = aes(x = Inf, y = Inf, label = round(mu$grp.mean, digits = 3)),
    hjust   = 1.05,
    vjust   = 1.5
  )+
  scale_color_brewer(palette="Dark2") +
  theme_bw()
dev.off()

### per region
library(reshape2)
beta_df$cpg <- rownames(beta_df)
betas_locations_m <- melt(data = beta_df[,c("cpg", "SUBJID", "Tissue", "Type")], id.vars = c('cpg','Type'),
                          variable.name = 'Variable', value.name = 'Proportion')
head(betas_locations_m)

library("dplyr")
mu <- betas_locations_m %>% 
  group_by(Variable , Type) %>%
  summarise(grp.mean = mean(Proportion))
mu

betas_locations_m$Type = factor(betas_locations_m$Type, levels=c('Promoter_Associated','Enhancer_Associated','Gene_Associated','Other'))
mu$Type = factor(mu$Type, levels=c('Promoter_Associated','Enhancer_Associated','Gene_Associated','Other'))

pdf(file = paste0('Plots/',"/VarPart_ann_main1.pdf"), w = 8, h = 3)
ggplot(betas_locations_m, aes(x = Proportion, fill=Type)) + 
  geom_density() +
  scale_fill_brewer(palette="Dark2") +
  geom_vline(data=mu, aes(xintercept=grp.mean, color=Type),
             linetype="dashed") +
  geom_text(
    size    = 5,
    data    = mu,
    mapping = aes(x = Inf, y = Inf, label = round(mu$grp.mean, digits = 3)),
    hjust   = 1.05,
    vjust   = 1.5
  )+
  scale_color_brewer(palette="Dark2") +
  facet_grid(Variable ~ Type)+
  theme_bw()
dev.off()

## plot together 
pdf(file = paste0('Plots/',"/VarPart_ann_main2.pdf"), w = 8, h = 3)
ggplot(betas_locations_m, aes(x = Proportion, fill=Variable)) + 
  geom_density() +
  scale_fill_brewer(palette="Dark2") +
  geom_vline(data=mu, aes(xintercept=grp.mean, color=Variable),
             linetype="dashed") +
  geom_text(
    size    = 5,
    data    = mu,
    mapping = aes(x = Inf, y = Inf, label = round(mu$grp.mean, digits = 3)),
    hjust   = 1.05,
    vjust   = 1.5
  )+
  scale_color_brewer(palette="Dark2") +
  facet_grid( ~ Type)+
  theme_bw()
dev.off()

pdf(file = paste0('Plots/',"/VarPart_ann_main2.v3.pdf"), w = 6, h = 3)
ggplot(betas_locations_m, aes(x = Proportion, fill=Type)) + 
  geom_density(alpha=0.7) +
  scale_fill_brewer(palette="Dark2") +
  geom_vline(data=mu, aes(xintercept=grp.mean, color=Type),
             linetype="dashed") +
  geom_text(
    size    = 5,
    data    = mu,
    mapping = aes(x = Inf, y = Inf, label = round(mu$grp.mean, digits = 3)),
    hjust   = 1.05,
    vjust   = 1.5
  )+
  scale_color_brewer(palette="Dark2") +
  facet_grid( ~ Variable)+
  theme_bw()
dev.off()

all <- betas_locations_m %>% 
  group_by(Variable) %>%
  summarise(grp.mean = mean(Proportion))
all

pdf(file = paste0('Plots/',"/VarPart_all_main2.pdf"), w = 4, h = 3)
ggplot(betas_locations_m, aes(x = Proportion, fill=Variable)) + 
  geom_density() +
  scale_fill_brewer(palette="Dark2") +
  geom_vline(data=all, aes(xintercept=grp.mean, color=Variable),
             linetype="dashed") +
  geom_text(
    size    = 5,
    data    = all,
    mapping = aes(x = Inf, y = Inf, label = round(all$grp.mean, digits = 3)),
    hjust   = 1.05,
    vjust   = 1.5
  )+
  scale_color_brewer(palette="Dark2") +
  theme_bw()
dev.off()

### enrichments #####
my_fisher <- function(type, betas_locations_m=betas_locations_m){
  #             DS      Not DS
  # type
  # other_types
  print(type)
  
  type_df <- betas_locations_m[betas_locations_m$Variable=='Tissue' & betas_locations_m$Proportion>0.5,]
  other_type <- betas_locations_m[betas_locations_m$Variable=='Tissue' & betas_locations_m$Proportion<=0.5,]
  type_diff <- nrow(type_df[type_df$Type == type,])
  type_notdiff <- nrow(type_df) - type_diff
  other_type_diff <- nrow(other_type[other_type$Type == type,])
  other_type_notdiff <- nrow(other_type) - other_type_diff
  
  ### test significance
  m <- matrix(c(type_diff, type_notdiff, other_type_diff, other_type_notdiff), 2,2, byrow = T)
  print(m)  
  
  m[is.na(m)] <- 0
  #m <- m[c(type,paste0('No ',type)),]
  rownames(m) <- c('Ind var', "Not Ind var")
  colnames(m) <- c(type,"Other")
  print(m)
  f <- fisher.test(m)
  print(f)
  return(list("f" = f, "m" = type_diff))
}
# Two-tailed Fisher test
#families <- as.vector(unique(chromhmm_cpgs$Lung$region_chromhmm))
families <- unique(betas_locations_m$Type)
fisher_results <- lapply(families, function(region) my_fisher(region))
names(fisher_results) <- families

read_data <- function(variables, data){ #Function to prepare data to plot and compute adjusted p value
  
  odds_ratio <- lapply(variables, function(type) data[[type]][['f']]$estimate)
  adj.P.Val <- p.adjust(sapply(variables, function(type) data[[type]][['f']]$p.value), method = "BH")
  CI_down <- lapply(variables, function(type) data[[type]][['f']]$conf.int[1])
  CI_up <- lapply(variables, function(type) data[[type]][['f']]$conf.int[2])
  sample_size <- lapply(variables, function(type) data[[type]][['m']])
  
  
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

results_fisher <- read_data(families, fisher_results)

results_fisher$sig[results_fisher$sig =="Sig"] <- "FDR < 0.05"
results_fisher$sig[results_fisher$sig =="not Sig"] <- "FDR >= 0.05"
results_fisher$sig <- factor(results_fisher$sig, levels=c("FDR >= 0.05", "FDR < 0.05"))
results_fisher$region <- factor(results_fisher$region, levels=rev(c("Other","Gene_Associated","Enhancer_Associated","Promoter_Associated")))
g <- ggplot(results_fisher, aes(x=log2(oddsRatio), y=region, colour=sig)) +
  geom_errorbar(aes(xmin=log2(CI_down), xmax=log2(CI_up)), width=.3) +
  geom_vline(xintercept = 0) +
  #xlim(0,20) + #Only for Lung to show the 0
  geom_point(size=3) + ylab('') + theme_bw() +
  #scale_colour_manual(values=c("#BFC0C0", "#E26D5C")) +
  scale_colour_manual(values=c("#E26D5C")) +
  xlab("log2(Odds ratio)") +
  theme(legend.title = element_blank(),
        axis.text.x = element_text(colour="black", size=13),
        axis.text.y = element_text(colour="black", size=14),
        legend.text = element_text(colour="black", size=13),
        axis.title.x = element_text(size=16),
        legend.spacing.y = unit(-0.05, "cm"),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.border = element_rect(colour = "black", linewidth=1)) +
  scale_y_discrete(breaks=c("Other","Gene_Associated","Enhancer_Associated","Promoter_Associated"),
                   labels=c("Intergenic","Gene-Associated","Enhancer-Associated","Promoter-Associated"))# + xlim(0, 3)
# pdf(file = paste0("~/marenostrum/Projects/GTEx_v8/Methylation/Plots/genomic_location_shared",'_',trait,".pdf"), w = 6, h = 3.5)
# print(g)
# dev.off()

#Plot sample sizes:

g2 <- ggplot(results_fisher) + geom_col(aes(sample_size, region), width = 0.6, col='#A39A92') +
  theme_classic() + xlab("Number of DMPs") + ylab("") +
  theme(legend.position = "none",
        axis.text.x = element_text(colour="black", size=13),
        axis.text.y=element_blank(),
        axis.title.x = element_text(size=16)) +
  scale_y_discrete(breaks=c("Other","Gene_Associated","Enhancer_Associated","Promoter_Associated"),
                   labels=c("Intergenic","Gene-Associated","Enhancer-Associated","Promoter-Associated")) +# + xlim(0, 3)
  scale_x_continuous(breaks=c(0, 70000, 200000)) #Only for lung
# pdf(file = paste0("~/marenostrum/Projects/GTEx_v8/Methylation/Plots/chromhmm/genomic_location_shared",'_',trait,"_sample_size.pdf"), w = 4, h = 3.5)
# print(g2)
# dev.off()

library(ggpubr)
p <- ggarrange(g, g2, labels = c("A", "B"),
               common.legend = TRUE, legend = "right", widths = c(0.8,0.3))
pdf(file = paste0(basepath, "/Projects/GTEx_v8/Methylation/Plots/Enrichment_location_tissue_variable.pdf"), w = 8, h = 3)
print(p)
dev.off()

# now we do the same plots but classified by CpG content 

mu_cpg$Type

pdf(file = paste0(basepath, '/Projects/GTEx_v8/Methylation/Plots/',"/VarPart_ann_main2.v3_cpg.pdf"), w = 6, h = 3)
ggplot(betas_locations_cpg, aes(x = Proportion, fill=Type)) + 
  geom_density(alpha=0.7) +
  scale_fill_brewer(palette="Dark2") +
  geom_vline(data=mu_cpg, aes(xintercept=grp.mean, color=Type),
             linetype="dashed") +
    geom_text(
      data = mu_cpg,
      aes(
        x = Inf,
        y = Inf,
        label = round(grp.mean, 3),
        color = Type, 
      ),
      hjust = 1.05,
      size = 5, 
      position = position_nudge()
    )+
  scale_color_brewer(palette="Dark2") +
  facet_grid( ~ Variable)+
  theme_bw()
dev.off()



library(reshape2)
annotation_cpg <- read.csv(paste0(basepath, "Projects/GTEx_v8/Methylation/Data/GPL21145_MethylationEPIC_15073387_v-1-0_processed.csv"))
beta_df$Type_cpg <- 'Open_sea'
beta_df$Type_cpg[rownames(beta_df) %in% annotation_cpg[annotation_cpg$Relation_to_UCSC_CpG_Island=="Island","IlmnID"]] <- "Island"
beta_df$Type_cpg[rownames(beta_df) %in% annotation_cpg[annotation_cpg$Relation_to_UCSC_CpG_Island %in% c("S_Shelf", "N_Shelf"),"IlmnID"]] <- "Shelf"
beta_df$Type_cpg[rownames(beta_df) %in% annotation_cpg[annotation_cpg$Relation_to_UCSC_CpG_Island %in% c("S_Shore", "N_Shore"),"IlmnID"]] <- "Shore"
table(beta_df$Type_cpg)
betas_locations_cpg <- melt(data = beta_df[,c("cpg", "SUBJID", "Tissue", "Type_cpg")], id.vars = c('cpg','Type_cpg'),
                          variable.name = 'Variable', value.name = 'Proportion')
colnames(betas_locations_cpg)[2] <- "Type"
head(betas_locations_cpg)

library("dplyr")
mu_cpg <- betas_locations_cpg %>% 
  group_by(Variable , Type) %>%
  summarise(grp.mean = mean(Proportion))
mu_cpg

betas_locations_cpg$Type = factor(betas_locations_cpg$Type, levels=c('Island','Shelf','Shore','Open_sea'))
mu_cpg$Type = factor(mu_cpg$Type, levels=c('Island','Shelf','Shore','Open_sea'))

families <- unique(betas_locations_cpg$Type)
fisher_results_cpg <- lapply(families, function(region) my_fisher(region, betas_locations_cpg))
names(fisher_results_cpg) <- families

results_fisher <- read_data(families, fisher_results_cpg)

results_fisher$sig[results_fisher$sig =="Sig"] <- "FDR < 0.05"
results_fisher$sig[results_fisher$sig =="not Sig"] <- "FDR >= 0.05"
results_fisher$sig <- factor(results_fisher$sig, levels=c("FDR >= 0.05", "FDR < 0.05"))
results_fisher$region <- factor(results_fisher$region, levels=rev(c('Island','Shelf','Shore','Open_sea')))
g <- ggplot(results_fisher, aes(x=log2(oddsRatio), y=region, colour=sig)) +
  geom_errorbar(aes(xmin=log2(CI_down), xmax=log2(CI_up)), width=.3) +
  geom_vline(xintercept = 0) +
  #xlim(0,20) + #Only for Lung to show the 0
  geom_point(size=3) + ylab('') + theme_bw() +
  #scale_colour_manual(values=c("#BFC0C0", "#E26D5C")) +
  scale_colour_manual(values=c("#E26D5C")) +
  xlab("log2(Odds ratio)") +
  theme(legend.title = element_blank(),
        axis.text.x = element_text(colour="black", size=13),
        axis.text.y = element_text(colour="black", size=14),
        legend.text = element_text(colour="black", size=13),
        axis.title.x = element_text(size=16),
        legend.spacing.y = unit(-0.05, "cm"),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.border = element_rect(colour = "black", linewidth=1)) #+
  scale_y_discrete(breaks=c("Other","Gene_Associated","Enhancer_Associated","Promoter_Associated"),
                   labels=c("Intergenic","Gene-Associated","Enhancer-Associated","Promoter-Associated"))# + xlim(0, 3)


  #Plot sample sizes:
  
  g2 <- ggplot(results_fisher) + geom_col(aes(sample_size, region), width = 0.6, col='#A39A92') +
    theme_classic() + xlab("Number of DMPs") + ylab("") +
    theme(legend.position = "none",
          axis.text.x = element_text(colour="black", size=13),
          axis.text.y=element_blank(),
          axis.title.x = element_text(size=16)) #+
    scale_y_discrete(breaks=c("Other","Gene_Associated","Enhancer_Associated","Promoter_Associated"),
                     labels=c("Intergenic","Gene-Associated","Enhancer-Associated","Promoter-Associated")) +# + xlim(0, 3)
    scale_x_continuous(breaks=c(0, 70000, 200000)) #Only for lung
  
  library(ggpubr)
  p <- ggarrange(g, g2, labels = c("A", "B"),
                 common.legend = TRUE, legend = "right", widths = c(0.8,0.3))
  pdf(file = paste0(basepath, "/Projects/GTEx_v8/Methylation/Plots/Enrichment_location_cpg_tissue_variable.pdf"), w = 8, h = 3)
  print(p)
  dev.off()


library(missMethyl)


GOenrichments <- list()
for (type in c('Promoter_Associated','Enhancer_Associated','Gene_Associated','Other')) {
  print(type)
  betas <- unique(betas_locations_m$cpg)
  GOenrichments[[type]] <- list()
    tryCatch(
      {res <- gometh(betas_locations_m$cpg[betas_locations_m$Type==type & betas_locations_m$Variable == 'SUBJID' & betas_locations_m$Proportion>=0.5], all.cpg=(betas_locations_m$cpg[betas_locations_m$Type==type]),
                     collection="GO", array.type="EPIC", sig.genes = TRUE)
      #res <- res[res$ONTOLOGY=="BP",]
      print(table(res$P.DE<0.05))
      if (sum(res$P.DE<0.05) > 0) {
        GOenrichments[[type]] <- res[res$P.DE<0.05,]
      }
      },  error=function(cond) {
        message("Error")
        message("Here's the original error message:")
        message(cond)
        # Choose a return value in case of error
        return(NA)
      })
}

GOenrichments <- do.call('rbind.data.frame', GOenrichments)
GOenrichments$location <- gsub('\\..*', '', rownames(GOenrichments))
GOenrichments$GO <- gsub('.*\\.', '', rownames(GOenrichments))
write.table(GOenrichments,'~/marenostrum/Projects/GTEx_v8/Methylation/Data/SupplTable1.1.Indiv_variable_GO_per_location.txt', sep = '\t', row.names = F,
            col.names = T, quote = F)

res <- gometh(betas_locations_m$cpg[betas_locations_m$Variable == 'SUBJID' & betas_locations_m$Proportion>=0.5], all.cpg=unique(betas_locations_m$cpg),
              collection="GO", array.type="EPIC", sig.genes = TRUE)
res$GO <- rownames(res)
write.table(res[res$P.DE<0.05,],'~/marenostrum/Projects/GTEx_v8/Methylation/Data/SupplTabl1_IndividualVariableGO.txt', sep = '\t', row.names = F,
            col.names = T, quote = F)

res <- res[res$ONTOLOGY=="BP",]

topgo <- topGSA(res, n=20)
topgo$Description <- rownames(topgo)

pdf(paste0("Plots/","GO_cpgs_tissue_variable.pdf"), width = 10, height = 5)
print(ggplot(data = topgo[topgo$FDR<0.1,], aes(x=DE/N, y = factor(TERM, levels=rev(TERM)),
                               color = -log10(FDR), size = DE)) +
        geom_point() + scale_color_gradient(low = "red", high = "blue") +
        theme_bw() +  ylab("") +  xlab("Gene Ratio") +
        ggtitle(" GO top all"))
dev.off()

### are individual variale, tissue shared? ######
sharing <- readRDS('~/marenostrum/Projects/GTEx_v8/Methylation/Data/Sharing_DMP.rds')
high_indiv_var_05 <- (betas_locations_m[betas_locations_m$Variable == 'SUBJID' & betas_locations_m$Proportion>0.5,])
table(high_indiv_var_05$Type)

head(sharing)
sharing_cpgs <- unique(sharing$CG[sharing$number>1 & sharing$trait=='EURv1'])
sum(high_indiv_var_05$cpg %in% sharing_cpgs)

#### Plot outliers 
# merge with gene info ##
Sys.time()
data_path <- "~/marenostrum/Projects/GTEx_v8/Methylation/Data/"
annotation <- read.delim(paste0(data_path, "Methylation_Epic_gene_promoter_enhancer_processed.txt"), sep = '\t', header = T)
Sys.time()

beta_df <- merge(beta_df, annotation, by.x='cpg', by.y='IlmnID', all.x=TRUE)

ggplot(beta_df, aes(x=Tissue, y=SUBJID)) +
  geom_bin2d(bins = 150) +
  scale_fill_continuous(type = "viridis") +
  theme_bw()#+
  #geom_text(aes(label=ifelse(SUBJID>0.8,as.character(UCSC_RefGene_Name),'')),hjust=0,vjust=0)

#### name genes in immune pathways #####
mhc <- read.table('Data/GO_MHC_Complex_assembly.txt', sep = '\t', row.names = NULL)
peptide <- read.table('Data/GO_antigen_processing.txt', sep = '\t', row.names = NULL)
beta_df$immune <- 'No'
beta_df$immune[beta_df$UCSC_RefGene_Name %in% c(mhc$V2, peptide$V2)] <- 'Immune'


library(ggrepel)
pdf(paste0("Plots/","Imune_cpgs_tissue_variable.dotplot.pdf"), width = 6, height = 5)
ggplot(beta_df[beta_df$SUBJID>0.7,], aes(x=Tissue, y=SUBJID)) +
  geom_point(alpha=0.6, colour = '#BFC0C0') +
  geom_point(data = beta_df[beta_df$SUBJID>0.7 & beta_df$immune == 'Immune',], aes(x=Tissue, y=SUBJID),alpha=0.8, colour='#058ED9') +
  geom_text_repel(aes(label=ifelse(immune=='Immune',as.character(UCSC_RefGene_Name),'')), max.overlaps = Inf) +
  theme_classic()
dev.off()

### plot only naming DMPs
sharing_dmps <- readRDS('~/marenostrum/Projects/GTEx_v8/Methylation/Data/Sharing_DMP.rds')
beta_df <- merge(beta_df, sharing_dmps[,c('CG','trait','number','dir')], by.x='cpg', by.y='CG', all.x=TRUE)

beta_df$label <- 'No'
beta_df$label[!is.na(beta_df$trait) & beta_df$SUBJID > 0.75] <- 'Yes'

colors_traits <- list('AGE'=c('#3D7CD0'),
                      'SEX2'=c('#3B734E'),
                      'EURv1'=c('#F0AE21'))

library(ggplot2)
library(ggExtra)
library(ggrepel)

for (trait in c('AGE','SEX2','EURv1')) {
  beta_df$label <- 'No'
  beta_df$label[beta_df$trait==trait & beta_df$SUBJID > 0.75] <- 'Yes'
  print(table( beta_df$label))
  pdf(paste0("Plots/","cpgs_ind_variable_",trait,".pdf"), width = 10, height = 5)
  p1 <- ggplot(beta_df[beta_df$trait==trait,], aes(x=Tissue, y=SUBJID, color=trait)) +
          geom_point() +
          scale_colour_manual(values=colors_traits[[trait]]) +
          theme_classic()+
          geom_text_repel(aes(label=ifelse(label=='Yes',as.character(UCSC_RefGene_Name),'')), max.overlaps = Inf, )
  p2 <- ggMarginal(p1, type="density")
  print(p2)
  dev.off()
}

genes_eur <- as.data.frame( table(beta_df[beta_df$trait==trait & beta_df$label=='Yes','UCSC_RefGene_Name']))

ggplot(genes_eur, aes(x=Freq, y=Var1)) +
  geom_bar(stat = 'identity', fill='dark grey') +
  theme_classic()+ xlab('Nº CpGs') + ylab('')

### compare expression ####

# Do enhancer contribute more to individual variance (age, sex and ancestry ) than promoters 

chuncks <- c(1:16)

#### reading varPart values ####
beta <- lapply(chuncks, function(chnk) readRDS(paste0('varPart/', chnk, "_chunck_var_part.rds")))
beta_df <- do.call("rbind",beta)



