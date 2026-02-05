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

to_plot_2 <- readRDS('~/marenostrum/Projects/GTEx_v8/Methylation/Data/correlations_to_plot_2.new.rds')
to_plot_1 <- readRDS('~/marenostrum/Projects/GTEx_v8/Methylation/Data/correlations_to_plot_1.new.rds')

#Barplot
### need to add color to facets
library(ggplot2)
to_plot_2 <- as.data.frame(to_plot_2[c(2:nrow(to_plot_2)),])
to_plot_2$V1 <- 100*as.numeric(to_plot_2$V1)
names(to_plot_2) <- c("N", "type", "Correlation", 'Tissue','Trait',"Number")
to_plot_2$type <- factor(to_plot_2$type, levels = c("Promoter", "Enhancer", "Gene Body"))

library(ggh4x)
g <- ggplot(to_plot_2[to_plot_2]) + geom_col(aes(type, N, fill=Correlation), width = 0.9) + xlab("") +
  # ylab("% of DEGs correlated with a DMP") +
  ylab("% of DEGs correlated with a DMP") +
  scale_fill_manual(values=c("#88CCEE", "#CC6677")) + theme_classic() +
  theme(axis.title.y = element_text(margin = margin(r = 2), size = 11),
        axis.text.y = element_text(size = 9, colour = "black"),
        axis.text.x = element_text(size = 9, colour = "black", angle = 90, vjust = 0.5)) +
  facet_nested(Trait ~ Tissue,scales = "free_y", independent = "y")
  #facet_grid2(vars(Trait), vars(Tissue), scales = "free_y", independent = "y")


# scale_fill_manual(values = rev(c(lighten("#88CCEE", amount = 0.5), "#88CCEE",
#                                  lighten("#CC6677", amount = 0.3), "#CC6677")),
#                   breaks = rev(levels(data$type)))

pdf("~/marenostrum/Projects/GTEx_v8/Methylation/Plots/DEGs_DMPs.New.pdf", width = 10, height = 4)
g
dev.off()



#% of DMPs correlated
to_plot_1 <- as.data.frame(to_plot_1[c(2:nrow(to_plot_1)),])
to_plot_1$V1 <- 100*as.numeric(to_plot_1$V1)
names(to_plot_1) <- c("N", "type", "Correlation", 'Tissue','Trait',"Number")
to_plot_1$type <- factor(to_plot_1$type, levels = c("Promoter", "Enhancer", "Gene Body"))

g <- ggplot(to_plot_1[to_plot_1$Trait != 'BMI',]) + geom_col(aes(type, N, fill=Correlation), width = 0.9) + xlab("") +
  # ylab("% of DMPs correlated with a DEG") +
  ylab("% of DMPs correlated with a DEG") +
  scale_fill_manual(values=c("#88CCEE", "#CC6677")) + theme_classic() +
  theme(axis.title.y = element_text(margin = margin(r = 2), size = 11),
        axis.text.y = element_text(size = 9, colour = "black"),
        axis.text.x = element_text(size = 9, colour = "black", angle = 90, vjust = 0.5)) +
  facet_nested(Trait ~ Tissue,scales = "free_y", independent = "y")

pdf("~/marenostrum/Projects/GTEx_v8/Methylation/Plots/DMPs_DEGs.New.pdf", width = 10, height = 4)
g
dev.off()

#### barplot total number correlated negative, positive
to_plot_1$Number <- as.numeric(to_plot_1$Number)
to_plot_1_s <- to_plot_1 %>% dplyr::group_by(Trait, Correlation) %>% dplyr::summarise(total=sum(Number))

traits_cols <- c('#C49122','#4B8C61','#70A0DF')
names(traits_cols) <- c('Ancestry','Sex','Age')
strip <- strip_themed(background_x = elem_list_rect(fill = traits_cols))

g <- ggplot(to_plot_1_s[to_plot_1_s$Trait != 'BMI',]) + geom_col(aes(Trait, total, fill=Correlation), width = 0.9) + xlab("") +
  # ylab("% of DMPs correlated with a DEG") +
  ylab("#DMPs correlated with a DEG") +
  scale_fill_manual(values=c("#88CCEE", "#CC6677")) + theme_classic() +
  theme(axis.title.y = element_text(margin = margin(r = 2), size = 11),
        axis.text.y = element_text(size = 9, colour = "black"),
        axis.text.x = element_text(size = 9, colour = "black")) #+
  #facet_wrap2(~ Trait, strip = strip)

pdf("~/marenostrum/Projects/GTEx_v8/Methylation/Plots/DMPs_DEGs.direction_all.pdf", width = 5, height = 4)
g
dev.off()


### heatmap number correlated ####
### read Correlations
get_corr <- function(tissue, trait){
  if(tissue %in% sex_tissues & trait == "SEX2"){
    NA
  }else{
    model <- readRDS(paste0(project_path, "Tissues/",tissue, "/",trait,'_Correlations_probes_genes_DEG_DMP.pnominal_ancestry_c.rds'))
    model[model$p.adj<0.05,]
  }
}
DMPs_cor <- lapply(c('EURv1','SEX2','AGE','BMI'), function(trait) lapply(tissues, function(tissue) get_corr(tissue, trait)))
names(DMPs_cor) <- c("Ancestry", "Sex", "Age", "BMI")
for(trait in c("Ancestry", "Sex", "Age", "BMI")){names(DMPs_cor[[trait]]) <- tissues}

get_pairs <- function(tissue, trait){
  if(tissue %in% sex_tissues & trait == "SEX2"){
    NA
  }else{
    model <- readRDS(paste0(project_path, "Tissues/",tissue, "/",trait,'_Correlations_probes_genes_DEG_DMP.pnominal_ancestry_c.rds'))
    model[!is.na(model$gene),]
  }
}
DMPs_DEGs <- lapply(c('EURv1','SEX2','AGE','BMI'), function(trait) lapply(tissues, function(tissue) get_pairs(tissue, trait)))
names(DMPs_DEGs) <- c("Ancestry", "Sex", "Age", "BMI")
for(trait in c("Ancestry", "Sex", "Age", "BMI")){names(DMPs_DEGs[[trait]]) <- tissues}


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

counts <- sapply(c("Ancestry", "Sex", "Age", "BMI"), function(trait)
  sapply(tissues, function(tissue) 
   nrow(DMPs_cor[[trait]][[tissue]])
  ))

counts <- sapply(c("Ancestry", "Sex", "Age", "BMI"), function(trait)
  sapply(tissues, function(tissue) 
    nrow(DMPs_DEGs[[trait]][[tissue]])
  ))

counts <- sapply(c("Ancestry", "Sex", "Age", "BMI"), function(trait)
  sapply(tissues, function(tissue) 
    nrow(DMPs_cor[[trait]][[tissue]])/nrow(DMPs_DEGs[[trait]][[tissue]])*100
  ))


# Row annotation --
names(n_samples) <- tissue_info$tissue_abbrv
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

rownames(counts) <- tissue_info$tissue_abbrv

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
without_NA <- replace(counts, is.na(counts), "")
counts_num <- matrix(as.numeric(counts),    # Convert to numeric matrix
       ncol = ncol(counts))
rownames(counts_num) <- tissue_info$tissue_abbrv
colnames(counts_num) <- c('Ancestry','Sex','Age','BMI')
ht <- Heatmap((counts_num),
              col= colorRamp2(seq(0,1000,length.out=9),
                              (brewer.pal(9, "BuPu"))),
              # col=colorRamp2( c(0,1,3000),
              #                 c("white","#F5F8FC","#1266b5") ),
              na_col = "white",
              cluster_rows = F,
              cluster_columns = F,
              # name = "DE signal",
              name = "#pairs probes",
              row_names_side = "left",
              column_names_side = "top",
              column_names_rot =  60,
              column_names_gp = gpar(fontsize = 12),
              column_names_max_height= unit(9, "cm"),
              row_names_gp = gpar(fontsize = 12),
              left_annotation = row_ha_left,
              cell_fun = function(j, i, x, y, width, height, fill) {
                grid.text(my_pretty_num_function(without_NA[i, j]), x, y, gp = gpar(fontsize = 12))}
              
)

head(counts)
counts[(counts)=='NULL'] <- NA

without_NA <- round(without_NA, digits = 0)
counts_num <- matrix(as.numeric(counts),    # Convert to numeric matrix
                     ncol = ncol(counts))
rownames(counts_num) <- tissue_info$tissue_abbrv
colnames(counts_num) <- c('Ancestry','Sex','Age','BMI')

without_NA <- round(counts_num, digits = 0)
without_NA <- replace(without_NA, is.na(without_NA), "")


ht <- Heatmap((counts_num),
              col= colorRamp2(seq(0,100,length.out=9),
                              (brewer.pal(9, "BuPu"))),
              # col=colorRamp2( c(0,1,3000),
              #                 c("white","#F5F8FC","#1266b5") ),
              na_col = "white",
              cluster_rows = F,
              cluster_columns = F,
              # name = "DE signal",
              name = "%DMPs associated to a DEG sig corr",
              row_names_side = "left",
              column_names_side = "top",
              column_names_rot =  60,
              column_names_gp = gpar(fontsize = 12),
              column_names_max_height= unit(9, "cm"),
              row_names_gp = gpar(fontsize = 12),
              left_annotation = row_ha_left,
              cell_fun = function(j, i, x, y, width, height, fill) {
                grid.text(my_pretty_num_function(without_NA[i, j]), x, y, gp = gpar(fontsize = 12))}
              
)


pdf("marenostrum/Projects/GTEx_v8/Methylation/Plots/Perc_corr_DMPs_DEGs_heatmap.pdf",
    width = 6, height = 4)
draw(ht)
    # heatmap_legend_side = "bottom")
dev.off()
