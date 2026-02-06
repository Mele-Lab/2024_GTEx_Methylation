#!/usr/bin/env Rscript
# @Author: Jose Miguel Ramirez; Adapted by Winona Oliveros
# @E-mail: winona.oliveros@bsc.es
# @Description: Code to Plor correlation results; general overview
# @software version: R=4.2.2

library(ComplexHeatmap)
library(RColorBrewer)
suppressPackageStartupMessages(library(circlize))

# ---- Data ----- ####
# Demographic traits ----
traits_cols <- c("Ancestry" = "#E69F00",
                 "Age" = "#56B4E9",
                 "Sex" =  "#009E73",
                 "BMI" = "#CC79A7")
traits <- names(traits_cols)

# Tissues ----
### read methylation results ####
#first_dir <- "~/marenostrum/"
first_dir <- "/home/mariasr/cluster/"
first_dir <- "/Users/mariasopenar/cluster/"
project_path <- paste0(first_dir, "Projects/GTEx_v8/Methylation/")

tissues <- c("Lung", "ColonTransverse", "Ovary", "Prostate", "BreastMammaryTissue", "MuscleSkeletal", "KidneyCortex", "Testis", "WholeBlood")
names <- c("Age", "Ancestry", "BMI", "Sex")

tissue_info <- readRDS(paste0(first_dir, "Projects/GTEx_v8/Methylation/Data/Tissue_info_whole.rds"))

tissue_info <- tissue_info[!grepl("BreastMammaryTissue_", tissue_info$tissue_ID),]
tissue_info <- tissue_info[tissue_info$tissue_ID %in% tissues,]

sex_tissues <- c('Ovary','Prostate','Testis')

tissues <- tissue_info$tissue_ID
# tissues cols
tissues_cols <- tissue_info$colcodes 
names(tissues_cols) <- tissue_info$tissue_abbrv

# Metadata ----
metadata <- lapply(tissues, function(tissue) {
  metadata_ind <- readRDS(paste0(project_path, "Tissues/",tissue, "/metadata.rds"))
  print("Reading Admixture results")
  admixture_ancestry <- read.table(paste0(first_dir,'/Projects/GTEx_v8/Methylation/admixture_inferred_ancestry.txt'))
  colnames(admixture_ancestry) <- c('SUBJID','AFRv1','EURv1','inferred_ancestry','AFRv2','EURv2')
  metadata_ind <- merge(metadata_ind, admixture_ancestry[,c("SUBJID","EURv1")], by='SUBJID')
  metadata_ind
})
names(metadata) <- tissues

for(tissue in sex_tissues[c(1)]){
  metadata[[tissue]]$SEX <- "2"
  metadata[[tissue]]$SEX <- as.factor(metadata[[tissue]]$SEX)
  metadata[[tissue]] <- metadata[[tissue]][, colnames(metadata$MuscleSkeletal)[c(1:8,10:21)]]
}
for(tissue in sex_tissues[c(2,3)]){
  metadata[[tissue]]$SEX <- "1"
  metadata[[tissue]]$SEX <- as.factor(metadata[[tissue]]$SEX)
  metadata[[tissue]] <- metadata[[tissue]][, colnames(metadata$MuscleSkeletal)[c(1:8,10:21)]]
}

to_plot_2 <- readRDS(paste0(project_path,'/Data/correlations_to_plot_2_ancestry_continous.new.rds')) # gene-based
to_plot_1 <- readRDS(paste0(project_path, "/Data/correlations_to_plot_1_ancestry_continous.new.rds")) # probe-based

#Barplot % of DMPs per trait and position
library(ggplot2)
to_plot_1 <- as.data.frame(to_plot_1[c(2:nrow(to_plot_1)),])
to_plot_1$V1 <- 100*as.numeric(to_plot_1$V1)
names(to_plot_1) <- c("N", "type", "Correlation", 'Tissue','Trait',"Number")
to_plot_1$type <- factor(to_plot_1$type, levels = c("Promoter", "Enhancer", "Gene Body"))

library(ggh4x)
library(dplyr)
library(ggplot2)

df_plot <- to_plot_1 %>%
  filter(Trait != "BMI") %>%
  mutate(Number = as.numeric(Number)) %>%
  group_by(Trait, type, Correlation) %>%
  summarise(Number = sum(Number), .groups = "drop")
trait_order <- c("Ancestry", "Sex", "Age")   # put your desired order here

df_plot <- df_plot %>%
  mutate(Trait = factor(Trait, levels = trait_order))

g <- ggplot(df_plot, aes(type, Number, fill = Correlation)) +
  geom_col(position = "fill", alpha = 0.8) +
  geom_text(aes(label = Number),
            position = position_fill(vjust = 0.5),
            size = 3) +
  xlab("") +
  ylab("% of DEGs correlated with a DMP\n in each direction") +
  scale_fill_manual(values = c("#88CCEE", "#CC6677")) +
  theme_classic() +
  theme(
    axis.title.y = element_text(margin = margin(r = 2), size = 11),
    axis.text.y  = element_text(size = 9, colour = "black"),
    axis.text.x  = element_text(size = 9, colour = "black", angle = 90, vjust = 0.5)
  ) +
  facet_wrap2(~Trait, strip = strip, nrow = 1)

g

pdf(paste0(project_path, "/Plots/DEGs_DMPs.New_ancestry_continous_padjusted.pdf"), width = 10, height = 4)
g
dev.off()


#Barplot % of DMPs per trait and position
library(ggplot2)
traits_cols <- c('#C49122','#4B8C61','#70A0DF','#A76595')
names(traits_cols) <- c("Ancestry", "Sex", "Age", "BMI")
strip <- strip_themed(background_x = elem_list_rect(fill = traits_cols[1:3]))
to_plot_1$Trait <- factor(to_plot_1$Trait, levels = c("Ancestry", "Sex", "Age", "BMI"))
colors_types <- c('#274c77','#6096ba','#a3cef1')
colors_types <- c('#320A28','#511730','#8E443D')
colors_types <- c('#223843','#C6D7C1','#86AC8F')
colors_types <- c('#a3cef1','#8E443D')

library(dplyr)

df_tissue <- to_plot_1 %>%
  filter(
    !is.na(Tissue),
    Tissue %in% c('Testis','Ovary','Prostate','Lung','ColonTransverse','BreastMammaryTissue'),
    Trait %in% c("Ancestry", "Sex", "Age")
  ) %>%
  mutate(
    N      = as.numeric(N),        # N is in % already in your plot object
    Number = as.numeric(Number),
    bg_row = ifelse(N > 0, 100 * Number / N, 0)  # reconstruct bg for that row
  ) %>%
  group_by(Trait, Tissue, Correlation) %>%
  summarise(
    num_corr = sum(Number, na.rm = TRUE),
    bg_total = sum(bg_row,  na.rm = TRUE),
    N = ifelse(bg_total > 0, 100 * num_corr / bg_total, 0),  # RECOMPUTED tissue-level %
    .groups = "drop"
  ) %>%
  group_by(Trait, Tissue) %>%
  mutate(
    Label_total = sum(num_corr, na.rm = TRUE),      # pos + neg counts per tissue+trait
    N_total = sum(N, na.rm = TRUE)                  # total % (pos+neg) for placing label
  ) %>%
  ungroup() %>%
  mutate(Trait = factor(Trait, levels = c("Ancestry", "Sex", "Age")))
library(ggplot2)
library(ggh4x)

g <- ggplot(df_tissue, aes(y = Tissue, x = N, fill = Correlation)) +
  geom_col(position = position_stack(reverse = TRUE), alpha = 0.8) +
  # one label per bar: use a separate data frame with 1 row per Tissue×Trait
  geom_text(
    aes(label = num_corr),
    position = position_stack(vjust = 0.5, reverse = TRUE),
    size = 3, size=3,hjust=1, angle=20
  ) +
  scale_x_continuous(breaks = seq(0, 100, 25), limits = c(0, 100)) +
  scale_fill_manual(values = c('#a3cef1','#8E443D')) +
  ylab('') + xlab('Proportion correlated DMPs (%)') +
  theme_classic() +
  facet_wrap2(~ Trait, strip = strip, nrow = 1)

g



pdf(paste0(project_path, "/Plots/DEGs_DMPs_proportion_tissue.New_ancestry_continous_padjusted.pdf"), width = 10, height = 4)
g
dev.off()




### heatmap number of DEGs correlated with DMPs ####
metadata <- lapply(tissues, function(tissue) {
  metadata_ind <- readRDS(paste0(project_path, "Tissues/",tissue, "/metadata.rds"))
  print("Reading Admixture results")
  admixture_ancestry <- read.table(paste0(project_path,'/admixture_inferred_ancestry.txt'))
  colnames(admixture_ancestry) <- c('SUBJID','AFRv1','EURv1','inferred_ancestry','AFRv2','EURv2')
  metadata_ind <- merge(metadata_ind, admixture_ancestry[,c("SUBJID","EURv1")], by='SUBJID')
  metadata_ind
})
names(metadata) <- tissues
n_samples <- sapply(tissues, function(tissue) nrow(metadata[[tissue]]))
names(n_samples) <- tissues



DEA_GTEx <- lapply(c(tissues), function(t) readRDS(paste0(project_path,'/Data/DEA/', t, "/", t, "DEA_results_Ancestry_continous_.results.rds")))
names(DEA_GTEx) <- tissues
DEA_GTEx$Lung$Ancestry_continous
## number of DEGs
meth_genes <- read.delim(paste0(project_path,'/Data/Methylation_Epic_gene_promoter_enhancer_processed.txt'))
get_pairs <- function(tissue, trait){
  if(tissue %in% sex_tissues & trait == "Sex"){
    NA
  }else{
    if(trait=="Ancestry"){
      trait_deg <- "Ancestry_continous"
    }else{
      trait_deg <- trait
    }
    length(DEA_GTEx[[tissue]][[trait_deg]]$gene.name.x[DEA_GTEx[[tissue]][[trait_deg]]$gene.name.x %in% meth_genes$UCSC_RefGene_Name & DEA_GTEx[[tissue]][[trait_deg]]$adj.P.Val<0.05])
  }
}

genes_DE_with_probe <- lapply(c("Ancestry", "Sex", "Age", "BMI"), function(trait) lapply(tissues, function(tissue) get_pairs(tissue, trait)))
names(genes_DE_with_probe) <- c("Ancestry", "Sex", "Age", "BMI")
for(trait in c("Ancestry", "Sex", "Age", "BMI")){names(genes_DE_with_probe[[trait]]) <- tissues}

get_corr <- function(tissue, trait){
  if(tissue %in% sex_tissues & trait == "SEX2"){
    NA
  }else{
    model <- readRDS(paste0(project_path, "Tissues/",tissue, "/",trait,'_Correlations_probes_genes_DEG_DMP.padjusted_ancestry_c.rds'))
    model[model$p.adj<0.05,]
  }
}
DMPs_cor <- lapply(c('EURv1','SEX2','AGE','BMI'), function(trait) lapply(tissues, function(tissue) get_corr(tissue, trait)))
names(DMPs_cor) <- c("Ancestry", "Sex", "Age", "BMI")
for(trait in c("Ancestry", "Sex", "Age", "BMI")){names(DMPs_cor[[trait]]) <- tissues}


counts<- sapply(c("Ancestry", "Sex", "Age", "BMI"), function(trait) {
  sapply(tissues, function(tissue) {
    
    # keep your sex tissue rule consistent with your earlier code
    if (tissue %in% sex_tissues && trait == "Sex") return(NA_real_)
    
    df_sig <- DMPs_cor[[trait]][[tissue]]           # already p.adj < 0.05
    den    <- genes_DE_with_probe[[trait]][[tissue]] # total DEGs with probes
    
    if (is.null(df_sig) || is.na(df_sig)[1] || nrow(df_sig) == 0) return(NA_real_)
    if (is.null(den)   || is.na(den)   || den == 0)               return(NA_real_)
    
    num <- length(unique(df_sig$gene))
    perc <- (num / den) * 100
    return(num)
  })
})

counts_p<- sapply(traits, function(trait) {
  sapply(tissues, function(tissue) {
    
    # keep your sex tissue rule consistent with your earlier code
    if (tissue %in% sex_tissues && trait == "Sex") return(NA_real_)
    
    df_sig <- DMPs_cor[[trait]][[tissue]]           # already p.adj < 0.05
    den    <- genes_DE_with_probe[[trait]][[tissue]] # total DEGs with probes
    
    if (is.null(df_sig) || is.na(df_sig)[1] || nrow(df_sig) == 0) return(NA_real_)
    if (is.null(den)   || is.na(den)   || den == 0)               return(NA_real_)
    
    num <- length(unique(df_sig$gene))
    perc <- (num / den) * 100
    return(perc)
  })
})


# Row annotation --
names(n_samples) <- tissue_info$tissue_abbrv
tissue_info <- tissue_info[tissues,]
row_ha_left <- HeatmapAnnotation("Samples" = anno_barplot(n_samples,
                                                          gp = gpar(fill = tissue_info$colcodes,
                                                                    col = tissue_info$colcodes),
                                                          border=F),
                                 gap = unit(0.25,"cm"),
                                 show_legend = T, 
                                 show_annotation_name = T,
                                 annotation_name_rot = 90,
                                 annotation_name_gp = gpar(fontsize = 10),
                                 which = "row")
# Column annotation --
traits_cols <- c('#C49122','#4B8C61','#70A0DF','#A76595')
names(traits_cols) <- c("Ancestry","Sex","Age","BMI")
column_ha_top <- HeatmapAnnotation("Traits" = as.character(1:4),
                                   col = list("Traits" = traits_cols),
                                   show_legend = T, show_annotation_name = F,
                                   simple_anno_size = unit(0.3,"cm"))

#rownames(counts) <- tissue_info$tissue_abbrv

# cell color % of tissue DEGs
my_pretty_num_function <- function(n){
  if(n==""){
    return(n)
  } else{
    prettyNum(n, big.mark = ",")
  }
}

head(counts)
counts[(counts)=='NULL'] <- NA
counts_num <- matrix(as.numeric(counts),    # Convert to numeric matrix
                     ncol = ncol(counts))
rownames(counts_num) <- tissues
colnames(counts_num) <- c('Ancestry','Sex','Age','BMI')
without_NA <- round(counts_num, digits = 2)
without_NA <- replace(without_NA, is.na(without_NA), "")

# counts_num <- matrix(as.numeric(counts),    # Convert to numeric matrix
#                      ncol = ncol(counts))
# rownames(counts_num) <- tissues
# colnames(counts_num) <- c('Ancestry','Sex','Age','BMI')
ht <- Heatmap((counts_p),
              col= colorRamp2(seq(0,30,length.out=9),
                              (brewer.pal(9, "BuPu"))),
              # col=colorRamp2( c(0,1,3000),
              #                 c("white","#F5F8FC","#1266b5") ),
              na_col = "white",
              cluster_rows = F,
              cluster_columns = F,
              # name = "DE signal",
              name = "% DEGs correlated",
              row_names_side = "left",
              column_names_side = "top",
              column_names_rot =  60,
              column_names_gp = gpar(fontsize = 12),
              column_names_max_height= unit(9, "cm"),
              row_names_gp = gpar(fontsize = 12),
              left_annotation = row_ha_left,
              cell_fun = function(j, i, x, y, width, height, fill) {
                grid.text(my_pretty_num_function((without_NA[i, j])), x, y, gp = gpar(fontsize = 12))}
              
)

pdf(paste0(project_path, "/Plots/Perc_explained_DEGs_heatmap.pajdusted_ancestry_continous.pdf"),
    width = 6, height = 4)
draw(ht)
# heatmap_legend_side = "bottom")
dev.off()










# #% of DMPs correlated
# to_plot_1 <- as.data.frame(to_plot_1[c(2:nrow(to_plot_1)),])
# to_plot_1$V1 <- 100*as.numeric(to_plot_1$V1)
# names(to_plot_1) <- c("N", "type", "Correlation", 'Tissue','Trait',"Number")
# to_plot_1$type <- factor(to_plot_1$type, levels = c("Promoter", "Enhancer", "Gene Body"))
# 
# g <- ggplot(to_plot_1[to_plot_1$Trait != 'BMI',]) + geom_col(aes(type, N, fill=Correlation), width = 0.9) + xlab("") +
#   # ylab("% of DMPs correlated with a DEG") +
#   ylab("% of DMPs correlated with a DEG") +
#   scale_fill_manual(values=c("#88CCEE", "#CC6677")) + theme_classic() +
#   theme(axis.title.y = element_text(margin = margin(r = 2), size = 11),
#         axis.text.y = element_text(size = 9, colour = "black"),
#         axis.text.x = element_text(size = 9, colour = "black", angle = 90, vjust = 0.5)) +
#   facet_nested(Trait ~ Tissue,scales = "free_y", independent = "y")
# 
# pdf("~/marenostrum/Projects/GTEx_v8/Methylation/Plots/DMPs_DEGs.New.pdf", width = 10, height = 4)
# g
# dev.off()
# 
# #### barplot total number correlated negative, positive
# to_plot_1$Number <- as.numeric(to_plot_1$Number)
# to_plot_1_s <- to_plot_1 %>% dplyr::group_by(Trait, Correlation) %>% dplyr::summarise(total=sum(Number))
# 
# traits_cols <- c('#C49122','#4B8C61','#70A0DF')
# names(traits_cols) <- c('Ancestry','Sex','Age')
# strip <- strip_themed(background_x = elem_list_rect(fill = traits_cols))
# 
# g <- ggplot(to_plot_1_s[to_plot_1_s$Trait != 'BMI',]) + geom_col(aes(Trait, total, fill=Correlation), width = 0.9) + xlab("") +
#   # ylab("% of DMPs correlated with a DEG") +
#   ylab("#DMPs correlated with a DEG") +
#   scale_fill_manual(values=c("#88CCEE", "#CC6677")) + theme_classic() +
#   theme(axis.title.y = element_text(margin = margin(r = 2), size = 11),
#         axis.text.y = element_text(size = 9, colour = "black"),
#         axis.text.x = element_text(size = 9, colour = "black")) +
#   facet_wrap2(~ Trait, strip = strip)
# 
# pdf("~/marenostrum/Projects/GTEx_v8/Methylation/Plots/DMPs_DEGs.direction_all.pdf", width = 5, height = 4)
# g
# dev.off()
# 
# 
# ### heatmap number correlated ####
# ### read Correlations
# get_corr <- function(tissue, trait){
#   if(tissue %in% sex_tissues & trait == "SEX2"){
#     NA
#   }else{
#     model <- readRDS(paste0(project_path, "Tissues/",tissue, "/",trait,'_Correlations_probes_genes_DEG_DMP.pnominal_ancestry_c.rds'))
#     model[model$p.adj<0.05,]
#   }
# }
# DMPs_cor <- lapply(c('EURv1','SEX2','AGE','BMI'), function(trait) lapply(tissues, function(tissue) get_corr(tissue, trait)))
# names(DMPs_cor) <- c("Ancestry", "Sex", "Age", "BMI")
# for(trait in c("Ancestry", "Sex", "Age", "BMI")){names(DMPs_cor[[trait]]) <- tissues}
# 
# get_pairs <- function(tissue, trait){
#   if(tissue %in% sex_tissues & trait == "SEX2"){
#     NA
#   }else{
#     model <- readRDS(paste0(project_path, "Tissues/",tissue, "/",trait,'_Correlations_probes_genes_DEG_DMP.pnominal_ancestry_c.rds'))
#     model[!is.na(model$gene),]
#   }
# }
# DMPs_DEGs <- lapply(c('EURv1','SEX2','AGE','BMI'), function(trait) lapply(tissues, function(tissue) get_pairs(tissue, trait)))
# names(DMPs_DEGs) <- c("Ancestry", "Sex", "Age", "BMI")
# for(trait in c("Ancestry", "Sex", "Age", "BMI")){names(DMPs_DEGs[[trait]]) <- tissues}
# 
# 
# metadata <- lapply(tissues, function(tissue) {
#   metadata_ind <- readRDS(paste0(project_path, "Tissues/",tissue, "/metadata.rds"))
#   print("Reading Admixture results")
#   admixture_ancestry <- read.table(paste0(project_path,'/admixture_inferred_ancestry.txt'))
#   colnames(admixture_ancestry) <- c('SUBJID','AFRv1','EURv1','inferred_ancestry','AFRv2','EURv2')
#   metadata_ind <- merge(metadata_ind, admixture_ancestry[,c("SUBJID","EURv1")], by='SUBJID')
#   metadata_ind
# })
# names(metadata) <- tissues
# n_samples <- sapply(tissues, function(tissue) nrow(metadata[[tissue]]))
# names(n_samples) <- tissues
# 
# counts <- sapply(c("Ancestry", "Sex", "Age", "BMI"), function(trait)
#   sapply(tissues, function(tissue) 
#    nrow(DMPs_cor[[trait]][[tissue]])
#   ))
# 
# counts <- sapply(c("Ancestry", "Sex", "Age", "BMI"), function(trait)
#   sapply(tissues, function(tissue) 
#     nrow(DMPs_DEGs[[trait]][[tissue]])
#   ))
# 
# counts <- sapply(c("Ancestry", "Sex", "Age", "BMI"), function(trait)
#   sapply(tissues, function(tissue) 
#     nrow(DMPs_cor[[trait]][[tissue]])/nrow(DMPs_DEGs[[trait]][[tissue]])*100
#   ))
# 
# 
# 
# 
# 
# 
# # Row annotation --
# names(n_samples) <- tissue_info$tissue_abbrv
# 
# row_ha_left <- HeatmapAnnotation("Samples" = anno_barplot(n_samples,
#                                                           gp = gpar(fill = tissue_info$colcodes,
#                                                                     col = tissue_info$colcodes),
#                                                           border=F),
#                                  gap = unit(0.25,"cm"),
#                                  show_legend = T, 
#                                  show_annotation_name = T,
#                                  annotation_name_rot = 90,
#                                  annotation_name_gp = gpar(fontsize = 10),
#                                  which = "row")
# # Column annotation --
# traits_cols <- c('#C49122','#4B8C61','#70A0DF','#A76595')
# names(traits_cols) <- c("Ancestry","Sex","Age","BMI")
# column_ha_top <- HeatmapAnnotation("Traits" = as.character(1:4),
#                                    col = list("Traits" = traits_cols),
#                                    show_legend = T, show_annotation_name = F,
#                                    simple_anno_size = unit(0.3,"cm"))
# 
# rownames(counts) <- tissue_info$tissue_abbrv
# 
# # cell color % of tissue DEGs
# my_pretty_num_function <- function(n){
#   if(n==""){
#     return(n)
#   } else{
#     prettyNum(n, big.mark = ",")
#   }
# }
# 
# head(counts)
# counts[(counts)=='NULL'] <- NA
# without_NA <- replace(counts, is.na(counts), "")
# counts_num <- matrix(as.numeric(counts),    # Convert to numeric matrix
#        ncol = ncol(counts))
# rownames(counts_num) <- tissue_info$tissue_abbrv
# colnames(counts_num) <- c('Ancestry','Sex','Age','BMI')
# ht <- Heatmap((counts_num),
#               col= colorRamp2(seq(0,1000,length.out=9),
#                               (brewer.pal(9, "BuPu"))),
#               # col=colorRamp2( c(0,1,3000),
#               #                 c("white","#F5F8FC","#1266b5") ),
#               na_col = "white",
#               cluster_rows = F,
#               cluster_columns = F,
#               # name = "DE signal",
#               name = "#pairs probes",
#               row_names_side = "left",
#               column_names_side = "top",
#               column_names_rot =  60,
#               column_names_gp = gpar(fontsize = 12),
#               column_names_max_height= unit(9, "cm"),
#               row_names_gp = gpar(fontsize = 12),
#               left_annotation = row_ha_left,
#               cell_fun = function(j, i, x, y, width, height, fill) {
#                 grid.text(my_pretty_num_function(without_NA[i, j]), x, y, gp = gpar(fontsize = 12))}
#               
# )
# 
# head(counts)
# counts[(counts)=='NULL'] <- NA
# 
# without_NA <- round(without_NA, digits = 0)
# counts_num <- matrix(as.numeric(counts),    # Convert to numeric matrix
#                      ncol = ncol(counts))
# rownames(counts_num) <- tissue_info$tissue_abbrv
# colnames(counts_num) <- c('Ancestry','Sex','Age','BMI')
# 
# without_NA <- round(counts_num, digits = 0)
# without_NA <- replace(without_NA, is.na(without_NA), "")
# 
# 
# ht <- Heatmap((counts_num),
#               col= colorRamp2(seq(0,100,length.out=9),
#                               (brewer.pal(9, "BuPu"))),
#               # col=colorRamp2( c(0,1,3000),
#               #                 c("white","#F5F8FC","#1266b5") ),
#               na_col = "white",
#               cluster_rows = F,
#               cluster_columns = F,
#               # name = "DE signal",
#               name = "%DMPs associated to a DEG sig corr",
#               row_names_side = "left",
#               column_names_side = "top",
#               column_names_rot =  60,
#               column_names_gp = gpar(fontsize = 12),
#               column_names_max_height= unit(9, "cm"),
#               row_names_gp = gpar(fontsize = 12),
#               left_annotation = row_ha_left,
#               cell_fun = function(j, i, x, y, width, height, fill) {
#                 grid.text(my_pretty_num_function(without_NA[i, j]), x, y, gp = gpar(fontsize = 12))}
#               
# )
# 
# 
# pdf("marenostrum/Projects/GTEx_v8/Methylation/Plots/Perc_corr_DMPs_DEGs_heatmap.pdf",
#     width = 6, height = 4)
# draw(ht)
#     # heatmap_legend_side = "bottom")
# dev.off()
# 
# 
# 
