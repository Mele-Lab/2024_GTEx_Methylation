#!/usr/bin/env Rscript
# @Author: Winona Oliveros Diez
# @E-mail: winn95@gmail.com
# @Description: Plot correlations, explore different threseholds of DMPs and other things
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
first_dir <- "~/marenostrum/"
project_path <- paste0(first_dir, "Projects/GTEx_v8/Methylation/")

tissues <- c("Lung", "ColonTransverse", "Ovary", "Prostate", "BreastMammaryTissue", "KidneyCortex", "Testis", "WholeBlood", "MuscleSkeletal")
names <- c("Age", "Ancestry", "BMI", "Sex")

tissue_info <- readRDS(paste0(first_dir, "Projects/GTEx_v8/Methylation/Data/Tissue_info_whole.rds"))

tissue_info <- tissue_info[!grepl("BreastMammaryTissue_", tissue_info$tissue_ID),]
tissue_info <- tissue_info[tissue_info$tissue_ID %in% tissues,]

sex_tissues <- c('Ovary','Prostate','Testis')

#tissues <- tissue_info$tissue_ID
# tissues cols
tissues_cols <- tissue_info$colcodes 
names(tissues_cols) <- tissue_info$tissue_abbrv

# Metadata ----
metadata <- lapply(tissues, function(tissue) {
  metadata_ind <- readRDS(paste0(project_path, "Tissues/",tissue, "/metadata.rds"))
  print("Reading Admixture results")
  admixture_ancestry <- read.table(paste0(project_path,'/admixture_inferred_ancestry.txt'))
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

#% of DMPs correlated
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
DMPs_co <- lapply(c('EURv1','SEX2','AGE','BMI'), function(trait) lapply(tissues, function(tissue) get_corr(tissue, trait)))
names(DMPs_Res) <- c("Ancestry", "Sex", "Age", "BMI")
for(trait in c("Ancestry", "Sex", "Age", "BMI")){names(DMPs_Res[[trait]]) <- tissues}


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

get_pairs <- function(tissue, trait){
  if(tissue %in% sex_tissues & trait == "SEX2"){
    NA
  }else{
    model <- readRDS(paste0(project_path, "Tissues/",tissue, "/",trait,'_Correlations_probes_genes_DEG_DMP.rds'))
    #model <- readRDS(paste0(project_path, "Tissues/",tissue, "/",trait,'_Correlations_probes_genes_DEG_DMP.rds'))
    #rownames(model[[trait]][model[[trait]]$adj.P.Val<0.05,])
    model[!is.na(model$gene),]
  }
}
DMPs_DEGs <- lapply(c('EURv1','SEX2','AGE','BMI'), function(trait) lapply(tissues, function(tissue) get_pairs(tissue, trait)))
names(DMPs_DEGs) <- c("Ancestry", "Sex", "Age", "BMI")
for(trait in c("Ancestry", "Sex", "Age", "BMI")){names(DMPs_DEGs[[trait]]) <- tissues}

DEA_GTEx <- lapply(c(tissues), function(t) readRDS(paste0(project_path,'/Data/DEA/', t, "/", t, "DEA_results_Ancestry_continous_.results.rds")))
names(DEA_GTEx) <- tissues
## number of DEGs
meth_genes <- read.delim(paste0(project_path,'/Data/Methylation_Epic_gene_promoter_enhancer_processed.txt'))
get_pairs <- function(tissue, trait){
  if(tissue %in% sex_tissues & trait == "Sex"){
    NA
  }else{
    length(DEA_GTEx[[tissue]][[trait]]$gene_name[DEA_GTEx[[tissue]][[trait]]$gene_name %in% meth_genes$UCSC_RefGene_Name & DEA_GTEx[[tissue]][[trait]]$adj.P.Val<0.05])
  }
}

genes_DE_with_probe <- lapply(c("Ancestry", "Sex", "Age", "BMI"), function(trait) lapply(tissues, function(tissue) get_pairs(tissue, trait)))
names(genes_DE_with_probe) <- c("Ancestry", "Sex", "Age", "BMI")
for(trait in c("Ancestry", "Sex", "Age", "BMI")){names(genes_DE_with_probe[[trait]]) <- tissues}

results_DML <- lapply(tissues, function(tis) 
  readRDS(paste0(paste0(project_path,"/Tissues/",tis,"/DML_results_5_PEERs_continous.rds"))))
names(results_DML) <- tissues

get_dmps_gene <- function(tissue, trait){
  if(tissue %in% sex_tissues & trait == "Sex"){
    NA
  }else{
    #sum(rownames(results_DML[[tissue]][[trait]][results_DML[[tissue]][[trait]]$adj.P.Val<0.05,]) %in% meth_genes$IlmnID)
    length(unique(meth_genes$UCSC_RefGene_Name[meth_genes$IlmnID %in% rownames(results_DML[[tissue]][[trait]][results_DML[[tissue]][[trait]]$adj.P.Val<0.05,])] ))
  }
}

DMPs_gene <- lapply(c('EURv1','SEX2','AGE','BMI'), function(trait) lapply(tissues, function(tissue) get_dmps_gene(tissue, trait)))
names(DMPs_gene) <- c("Ancestry", "Sex", "Age", "BMI")
for(trait in c("Ancestry", "Sex", "Age", "BMI")){names(DMPs_gene[[trait]]) <- tissues}


n_samples <- sapply(tissues, function(tissue) nrow(metadata[[tissue]]))
names(n_samples) <- tissues

counts <- sapply(c("Ancestry", "Sex", "Age", "BMI"), function(trait)
  sapply(tissues, function(tissue) 
    if (tissue %in% sex_tissues & trait == "Sex") {
      NA
    }else if (!nrow(DMPs_DEGs[[trait]][[tissue]])>0) {
      NA
    } else {
    length(unique(DMPs_DEGs[[trait]][[tissue]]$gene))/DMPs_gene[[trait]][[tissue]]}
  ))

counts <- sapply(c("Ancestry", "Sex", "Age", "BMI"), function(trait)
  sapply(tissues, function(tissue) 
    nrow(DMPs_Res[[trait]][[tissue]])
  ))

counts <- sapply(c("Ancestry", "Sex", "Age", "BMI"), function(trait)
  sapply(tissues, function(tissue) 
    nrow(DMPs_DEGs[[trait]][[tissue]])
  ))

counts_p <- sapply(c("Ancestry", "Sex", "Age", "BMI"), function(trait)
  sapply(tissues, function(tissue) 
    if (tissue %in% sex_tissues & trait == "Sex") {
      NA
    }else if (!nrow(DMPs_DEGs[[trait]][[tissue]])>0) {
      NA
    } else {
      length(unique(DMPs_DEGs[[trait]][[tissue]]$gene))/(genes_DE_with_probe[[trait]][[tissue]])*100
    }
  ))

counts <- sapply(c("Ancestry", "Sex", "Age", "BMI"), function(trait)
  sapply(tissues, function(tissue) 
    if (tissue %in% sex_tissues & trait == "Sex") {
      NA
    }else if (!nrow(DMPs_DEGs[[trait]][[tissue]])>0) {
      NA
    } else {
      length(unique(DMPs_DEGs[[trait]][[tissue]]$gene))
    }
  ))

counts_p <- sapply(c("Ancestry", "Sex", "Age", "BMI"), function(trait)
  sapply(tissues, function(tissue) 
    if (tissue %in% sex_tissues & trait == "Sex") {
      NA
    }else if (!nrow(signif[[trait]][[tissue]])>0) {
      NA
    } else {
      length(unique(signif[[trait]][[tissue]]$gene))/(genes_DE_with_probe[[trait]][[tissue]])*100
    }
  ))

counts <- sapply(c("Ancestry", "Sex", "Age", "BMI"), function(trait)
  sapply(tissues, function(tissue) 
    if (tissue %in% sex_tissues & trait == "Sex") {
      NA
    }else if (!nrow(signif[[trait]][[tissue]])>0) {
      NA
    } else {
      length(unique(signif[[trait]][[tissue]]$gene))
    }
  ))




# Row annotation --
#names(n_samples) <- tissue_info$tissue_abbrv
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
              col= colorRamp2(seq(0,50,length.out=9),
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

pdf("~/marenostrum/Projects/GTEx_v8/Methylation/Plots/Perc_explained_DEGs_heatmap.pnominal.pdf",
    width = 6, height = 4)
draw(ht)
# heatmap_legend_side = "bottom")
dev.off()

pdf("~/marenostrum/Projects/GTEx_v8/Methylation/Plots/Perc_correlated_DEGs_heatmap.FDR.number.pdf",
    width = 6, height = 4)
draw(ht)
# heatmap_legend_side = "bottom")
dev.off()

pdf("~/marenostrum/Projects/GTEx_v8/Methylation/Plots/Perc_correlated_DEGs_heatmap.pnominal.number.pdf",
    width = 6, height = 4)
draw(ht)
# heatmap_legend_side = "bottom")
dev.off()

pdf("~/marenostrum/Projects/GTEx_v8/Methylation/Plots/Number_DEGs_pairs.pnominal.pdf",
    width = 6, height = 4)
draw(ht)
# heatmap_legend_side = "bottom")
dev.off()

head(counts)
counts[(counts)=='NULL'] <- NA

counts_num <- matrix(as.numeric(counts),    # Convert to numeric matrix
                     ncol = ncol(counts))
rownames(counts_num) <- tissues
colnames(counts_num) <- c('Ancestry','Sex','Age','BMI')

without_NA <- round(counts_num, digits = 0)
without_NA <- replace(without_NA, is.na(without_NA), "")


ht <- Heatmap((counts_num),
              col= colorRamp2(seq(0,240,length.out=9),
                              (brewer.pal(9, "BuPu"))),
              # col=colorRamp2( c(0,1,3000),
              #                 c("white","#F5F8FC","#1266b5") ),
              na_col = "white",
              cluster_rows = F,
              cluster_columns = F,
              # name = "DE signal",
              name = "# DEG explained by Methylation",
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


pdf("~/marenostrum/Projects/GTEx_v8/Methylation/Plots/Number_explained_DEGs_heatmap.pnominal.pdf",
    width = 6, height = 4)
draw(ht)
# heatmap_legend_side = "bottom")
dev.off()

### number of DMPs per gene
get_corr <- function(tissue, trait){
  if(tissue %in% sex_tissues & trait == "SEX2"){
    NA
  }else{
    model <- readRDS(paste0(project_path, "Tissues/",tissue, "/",trait,'_Correlations_probes_genes_DEG_DMP.rds'))
    #rownames(model[[trait]][model[[trait]]$adj.P.Val<0.05,])
    model
  }
}
DMPs_Res_fdr <- lapply(c('EURv1','SEX2','AGE','BMI'), function(trait) lapply(tissues, function(tissue) get_corr(tissue, trait)))
names(DMPs_Res_fdr) <- c("Ancestry", "Sex", "Age", "BMI")
for(trait in c("Ancestry", "Sex", "Age", "BMI")){names(DMPs_Res_fdr[[trait]]) <- tissues}


get_signif <- function(tissue, trait){
  if(tissue %in% sex_tissues & trait == "Sex"){
    NA
  }else{
    DMPs_Res_fdr[[trait]][[tissue]][DMPs_Res_fdr[[trait]][[tissue]]$p.adj<0.05,]
  }
}
signif_fdr <- lapply(c("Ancestry", "Sex", "Age", "BMI"), function(trait) lapply(tissues, function(tissue) get_signif(tissue, trait)))
names(signif_fdr) <- c("Ancestry", "Sex", "Age", "BMI")
for(trait in c("Ancestry", "Sex", "Age", "BMI")){names(signif_fdr[[trait]]) <- tissues}

get_corr <- function(tissue, trait){
  if(tissue %in% sex_tissues & trait == "SEX2"){
    NA
  }else{
    model <- readRDS(paste0(project_path, "Tissues/",tissue, "/",trait,'_Correlations_probes_genes_DEG_DMP.pnominal.rds'))
    #rownames(model[[trait]][model[[trait]]$adj.P.Val<0.05,])
    model
  }
}
DMPs_Res <- lapply(c('EURv1','SEX2','AGE','BMI'), function(trait) lapply(tissues, function(tissue) get_corr(tissue, trait)))
names(DMPs_Res) <- c("Ancestry", "Sex", "Age", "BMI")
for(trait in c("Ancestry", "Sex", "Age", "BMI")){names(DMPs_Res[[trait]]) <- tissues}


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

all_cor <- unlist(signif, recursive = FALSE)
all_cor_df <- do.call("rbind", all_cor)
all_cor_df[,c('Trait','Tissue','Gene')] <- stringr::str_split_fixed(rownames(all_cor_df), '\\.', 3)
all_cor_df$label <- paste0(all_cor_df$gene,':',all_cor_df$probe, ':', all_cor_df$Trait, ':',all_cor_df$Tissue)

write.table(all_cor_df[,c("gene","probe","cor","p.val","p.adj","class","Trait","Tissue")], '~/marenostrum/Projects/GTEx_v8/Methylation/Data/Correlations_all_traits_tissues.pnominal.txt', sep = '\t', 
            col.names = T, row.names = F)


all_cor_fdr <- unlist(signif_fdr, recursive = FALSE)
all_cor_fdr <- do.call("rbind", all_cor_fdr)
all_cor_fdr[,c('Trait','Tissue','Gene')] <- stringr::str_split_fixed(rownames(all_cor_fdr), '\\.', 3)
all_cor_fdr$label <- paste0(all_cor_fdr$gene,':',all_cor_fdr$probe, ':', all_cor_fdr$Trait, ':',all_cor_fdr$Tissue)

res_list <- list(fdr=all_cor_fdr$label, pnominal=all_cor_df$label)
m1 = make_comb_mat(res_list)
cs = comb_size(m1)
upset <- UpSet(m1,
               top_annotation = HeatmapAnnotation(
  "Intersections" = anno_barplot(cs, 
                                       ylim = c(0, max(cs)*1.1),
                                       border = FALSE, 
                                       gp = gpar(fill = "black"), 
                                       height = unit(4, "cm")
  ), 
  annotation_name_side = "left", 
  annotation_name_rot = 90),)
upset <- draw(upset)
od = column_order(upset)
decorate_annotation("Intersections", {
  grid.text(cs[od], x = seq_along(cs), y = unit(cs[od], "native") + unit(2, "pt"), 
            default.units = "native", just = c("left", "bottom"), 
            gp = gpar(fontsize = 10, col = "#404040"), rot = 45)
})

pdf("~/marenostrum/Projects/GTEx_v8/Methylation/Plots/comparison_pnominal_fdr_upset.pdf",
    width = 6, height = 4)
upset <- draw(upset)
od = column_order(upset)
decorate_annotation("Intersections", {
  grid.text(cs[od], x = seq_along(cs), y = unit(cs[od], "native") + unit(2, "pt"), 
            default.units = "native", just = c("left", "bottom"), 
            gp = gpar(fontsize = 10, col = "#404040"), rot = 45)
})
# heatmap_legend_side = "bottom")
dev.off()

### number dmps per gene 
all_cor_df$type <- 'Nominal'
all_cor_df$type[all_cor_df$label %in% all_cor_fdr$label] <- 'FDR'

library(dplyr)
grouped <- all_cor_df %>% group_by(type,Tissue,gene,Trait) %>% summarise(n= n())

par(
  mfrow=c(1,2),
  mar=c(4,4,1,0)
)
pdf("~/marenostrum/Projects/GTEx_v8/Methylation/Plots/Number_DMPS_per_DEG_nominal_fdr.pdf",
    width = 6, height = 4)
par(
  mfrow=c(1,2),
  mar=c(4,4,1,0)
)
hist(grouped$n[grouped$type=='FDR'], breaks=40 , xlim=c(0,40) , col=rgb(1,0,0,0.5) , xlab="# DMPs/DEG FDR" , ylab="# DEGs" , main="" )
hist(grouped$n[grouped$type=='Nominal'], breaks=40 , xlim=c(0,40) , col=rgb(0,0,1,0.5) , xlab="# DMPs/DEG Nominal" , ylab="" , main="")
dev.off()

pdf("~/marenostrum/Projects/GTEx_v8/Methylation/Plots/Number_DMPS_per_DEG_nominal_fdr.barplot.pdf",
    width = 5, height = 3.5)
ggplot(grouped, aes(x=n)) + 
  geom_bar(stat = "count", width=0.5) + facet_wrap(~type) + ylab('# DEG') + xlab('# DMPs/DEG FDR')
dev.off()

### direction 
all_cor_fdr_tested <- unlist(DMPs_Res_fdr, recursive = FALSE)
all_cor_fdr_tested <- do.call("rbind", all_cor_fdr_tested)
all_cor_fdr_tested[,c('Trait','Tissue','Gene')] <- stringr::str_split_fixed(rownames(all_cor_fdr_tested), '\\.', 3)
all_cor_fdr_tested$label <- paste0(all_cor_fdr_tested$gene,':',all_cor_fdr_tested$probe, ':', all_cor_fdr_tested$Trait, ':',all_cor_fdr_tested$Tissue)

to_plot_1 <- c(1:6)
to_plot_2 <- c(1:6)
for (tissue in tissues) {
  for (trait in c("Ancestry", "Sex", "Age", "BMI")) {
    if(tissue %in% sex_tissues & trait == "Sex"){
      next
    } else {
      if(option==2){
        next
      } else {
        # all
        # pro_all <- DMPs_Res[[trait]][[tissue]][!is.na(DMPs_Res[[trait]][[tissue]]$gene),]
        # #pro_deg <- pro_all[pro_all$gene %in% degs,] #DEG-probe pairs
        # bg <- length(unique(pro_all$gene)) # DEGs that have at least one probe in promoters
        # corr <- pro_all[pro_all$p.adj<0.05,]
        # 
        # ### gene numbers 
        # # to_plot_2 <- rbind(to_plot_2,c(length(unique(corr[corr$cor>0,"gene"]))/bg, "Positive", tissue, trait,length(unique(corr[corr$cor>0,"gene"]))),
        # #                    c(length(unique(corr[corr$cor<0,"gene"]))/bg, "Negative", tissue, trait,length(unique(corr[corr$cor<0,"gene"]))))
        # 
        # bg <- length(unique(pro_all$probe)) # DMPs that are associated to a gene
        # ### probe numbers 
        # to_plot_1 <- rbind(to_plot_1,c(length(unique(corr[corr$cor>0,"probe"]))/bg, "Positive", tissue, trait,length(unique(corr[corr$cor>0,"probe"]))),
        #                    c(length(unique(corr[corr$cor<0,"probe"]))/bg, "Negative", tissue, trait,length(unique(corr[corr$cor<0,"probe"]))))
        # 
        #Promoters
        pro_all_fdr <- DMPs_Res_fdr[[trait]][[tissue]][!is.na(DMPs_Res_fdr[[trait]][[tissue]]$gene),]
        if (nrow(pro_all_fdr)>0) {
        pro_all_fdr$label <- paste0(pro_all_fdr$gene,':',pro_all_fdr$probe, ':', trait, ':',tissue)
        fdr <- pro_all_fdr$label
        } else {
          fdr <- c(NA)
        }
        pro_all <- DMPs_Res[[trait]][[tissue]][!is.na(DMPs_Res[[trait]][[tissue]]$gene) & DMPs_Res[[trait]][[tissue]]$class=="promoter",]
        if (nrow(pro_all)>0) {
          pro_all$label <- paste0(pro_all$gene,':',pro_all$probe, ':', trait, ':',tissue)
          #pro_deg <- pro_all[pro_all$gene %in% degs,] #DEG-probe pairs
          bg <- length(unique(pro_all$probe[pro_all$label %in% fdr])) 
          corr <- pro_all[pro_all$p.adj<0.05 & pro_all$label %in% fdr,]
          
          ### fdr numbers
          to_plot_2 <- rbind(to_plot_2,c(length(unique(corr[corr$cor>0,"probe"]))/bg, "Promoter", "Positive", tissue, trait,length(unique(corr[corr$cor>0,"probe"]))),
                             c(length(unique(corr[corr$cor<0,"probe"]))/bg, "Promoter", "Negative", tissue, trait,length(unique(corr[corr$cor<0,"probe"]))))
          
          bg <- length(unique(pro_all$probe[!pro_all$label %in% fdr])) 
          corr <- pro_all[!pro_all$label %in% fdr & pro_all$p.adj<0.05,]
          ### nominal numbers
          to_plot_1 <- rbind(to_plot_1,c(length(unique(corr[corr$cor>0,"probe"]))/bg, "Promoter", "Positive", tissue, trait,length(unique(corr[corr$cor>0,"probe"]))),
                             c(length(unique(corr[corr$cor<0,"probe"]))/bg, "Promoter", "Negative", tissue, trait,length(unique(corr[corr$cor<0,"probe"]))))
          
        } 
      
        #Enhancers
        pro_all <- DMPs_Res[[trait]][[tissue]][!is.na(DMPs_Res[[trait]][[tissue]]$gene) & DMPs_Res[[trait]][[tissue]]$class=="enhancer",]
        #pro_deg <- pro_all[pro_all$gene %in% degs,] #DEG-probe pairs
        if (nrow(pro_all)>0) {
          pro_all$label <- paste0(pro_all$gene,':',pro_all$probe, ':', trait, ':',tissue)
          bg <- length(unique(pro_all$probe[pro_all$label %in% fdr])) 
          corr <- pro_all[pro_all$p.adj<0.05 & pro_all$label %in% fdr,]
  
          ### gene numbers
          to_plot_2 <- rbind(to_plot_2,c(length(unique(corr[corr$cor>0,"probe"]))/bg, "Enhancer", "Positive", tissue, trait,length(unique(corr[corr$cor>0,"probe"]))),
                             c(length(unique(corr[corr$cor<0,"probe"]))/bg, "Enhancer", "Negative", tissue, trait,length(unique(corr[corr$cor<0,"probe"]))))
  
          bg <- length(unique(pro_all$probe[!pro_all$label %in% fdr])) 
          corr <- pro_all[!pro_all$label %in% fdr & pro_all$p.adj<0.05,]
          ### probe numbers
          to_plot_1 <- rbind(to_plot_1,c(length(unique(corr[corr$cor>0,"probe"]))/bg, "Enhancer", "Positive", tissue, trait,length(unique(corr[corr$cor>0,"probe"]))),
                             c(length(unique(corr[corr$cor<0,"probe"]))/bg, "Enhancer", "Negative", tissue, trait,length(unique(corr[corr$cor<0,"probe"]))))

        }
        #Gene Body
        pro_all <- DMPs_Res[[trait]][[tissue]][!is.na(DMPs_Res[[trait]][[tissue]]$gene) & DMPs_Res[[trait]][[tissue]]$class=="gene_body",]
        #pro_deg <- pro_all[pro_all$gene %in% degs,] #DEG-probe pairs
        if (nrow(pro_all)>0) {
          pro_all$label <- paste0(pro_all$gene,':',pro_all$probe, ':', trait, ':',tissue)
          bg <- length(unique(pro_all$probe[pro_all$label %in% fdr])) 
          corr <- pro_all[pro_all$p.adj<0.05 & pro_all$label %in% fdr,]
  
          ### gene numbers
          to_plot_2 <- rbind(to_plot_2,c(length(unique(corr[corr$cor>0,"probe"]))/bg, "Gene Body", "Positive", tissue, trait,length(unique(corr[corr$cor>0,"probe"]))),
                             c(length(unique(corr[corr$cor<0,"probe"]))/bg, "Gene Body", "Negative", tissue, trait,length(unique(corr[corr$cor<0,"probe"]))))
  
          
          bg <- length(unique(pro_all$probe[!pro_all$label %in% fdr])) 
          corr <- pro_all[!pro_all$label %in% fdr & pro_all$p.adj<0.05,]
          ### probe numbers
          to_plot_1 <- rbind(to_plot_1,c(length(unique(corr[corr$cor>0,"probe"]))/bg, "Gene Body", "Positive", tissue, trait,length(unique(corr[corr$cor>0,"probe"]))),
                             c(length(unique(corr[corr$cor<0,"probe"]))/bg, "Gene Body", "Negative", tissue, trait,length(unique(corr[corr$cor<0,"probe"]))))

        }
      }
    }
  }
}



to_plot_1 <- as.data.frame(to_plot_1[c(2:nrow(to_plot_1)),])
to_plot_1$V1 <- 100*as.numeric(to_plot_1$V1)
names(to_plot_1) <- c("N", "type", "Correlation", 'Tissue','Trait',"Number")
to_plot_1$type <- factor(to_plot_1$type, levels = c("Promoter", "Enhancer", "Gene Body"))

#pnominal
library(ggh4x)
g <- ggplot(to_plot_1[to_plot_1$Trait != 'BMI',]) + geom_col(aes(type, N, fill=Correlation), width = 0.9) + xlab("") +
  # ylab("% of DMPs correlated with a DEG") +
  ylab("% of DMPs correlated with a DEG") +
  scale_fill_manual(values=c("#88CCEE", "#CC6677")) + theme_classic() +
  theme(axis.title.y = element_text(margin = margin(r = 2), size = 11),
        axis.text.y = element_text(size = 9, colour = "black"),
        axis.text.x = element_text(size = 9, colour = "black", angle = 90, vjust = 0.5)) +
  facet_nested(Trait ~ Tissue,scales = "free_y", independent = "y")

pdf("~/marenostrum/Projects/GTEx_v8/Methylation/Plots/DMPs_DEGs.nominal_direction.pdf", width = 10, height = 4)
g
dev.off()

to_plot_2 <- as.data.frame(to_plot_2[c(2:nrow(to_plot_2)),])
to_plot_2$V1 <- 100*as.numeric(to_plot_2$V1)
names(to_plot_2) <- c("N", "type", "Correlation", 'Tissue','Trait',"Number")
to_plot_2$type <- factor(to_plot_2$type, levels = c("Promoter", "Enhancer", "Gene Body"))

#FDR
library(ggh4x)
g <- ggplot(to_plot_2[to_plot_2$Trait != 'BMI',]) + geom_col(aes(type, N, fill=Correlation), width = 0.9) + xlab("") +
  # ylab("% of DMPs correlated with a DEG") +
  ylab("% of DMPs correlated with a DEG") +
  scale_fill_manual(values=c("#88CCEE", "#CC6677")) + theme_classic() +
  theme(axis.title.y = element_text(margin = margin(r = 2), size = 11),
        axis.text.y = element_text(size = 9, colour = "black"),
        axis.text.x = element_text(size = 9, colour = "black", angle = 90, vjust = 0.5)) +
  facet_nested(Trait ~ Tissue,scales = "free_y", independent = "y")

pdf("~/marenostrum/Projects/GTEx_v8/Methylation/Plots/DMPs_DEGs.FDR_direction.pdf", width = 10, height = 4)
g
dev.off()

### plot merging all tissues 
to_plot_1 <- c(1:5)
to_plot_2 <- c(1:5)
  for (trait in c("Ancestry", "Sex", "Age", "BMI")) {

        #Promoters
        pro_all_fdr <- do.call('rbind.data.frame',DMPs_Res_fdr[[trait]])
        pro_all_fdr <- pro_all_fdr[!is.na(pro_all_fdr$gene),]
        pro_all_fdr$tissue <- gsub('\\..*','',rownames(pro_all_fdr))
        if (nrow(pro_all_fdr)>0) {
          pro_all_fdr$label <- paste0(pro_all_fdr$gene,':',pro_all_fdr$probe, ':', trait, ':',pro_all_fdr$tissue)
          fdr <- pro_all_fdr$label
        } else {
          fdr <- c(NA)
        }
        pro_all <- do.call('rbind.data.frame',DMPs_Res[[trait]])
        pro_all <- pro_all[!is.na(pro_all$gene) & pro_all$class=="promoter",]
        pro_all$tissue <- gsub('\\..*','',rownames(pro_all))
        if (nrow(pro_all)>0) {
          pro_all$label <- paste0(pro_all$gene,':',pro_all$probe, ':', trait, ':',pro_all$tissue)
          #pro_deg <- pro_all[pro_all$gene %in% degs,] #DEG-probe pairs
          bg <- length(unique(pro_all$label[pro_all$label %in% fdr])) 
          corr <- pro_all[pro_all$p.adj<0.05 & pro_all$label %in% fdr,]
          
          ### fdr numbers
          to_plot_2 <- rbind(to_plot_2,c(length(unique(corr[corr$cor>0,"label"]))/bg, "Promoter", "Positive", trait,length(unique(corr[corr$cor>0,"label"]))),
                             c(length(unique(corr[corr$cor<0,"label"]))/bg, "Promoter", "Negative", trait,length(unique(corr[corr$cor<0,"label"]))))
          
          bg <- length(unique(pro_all$label))#[!pro_all$label %in% fdr])) 
          corr <- pro_all[pro_all$p.adj<0.05,] # !pro_all$label %in% fdr & 
          ### nominal numbers
          to_plot_1 <- rbind(to_plot_1,c(length(unique(corr[corr$cor>0,"label"]))/bg, "Promoter", "Positive", trait,length(unique(corr[corr$cor>0,"label"]))),
                             c(length(unique(corr[corr$cor<0,"label"]))/bg, "Promoter", "Negative", trait,length(unique(corr[corr$cor<0,"label"]))))
          
        } 
        
        #Enhancers
        pro_all <- do.call('rbind.data.frame',DMPs_Res[[trait]])
        pro_all <- pro_all[!is.na(pro_all$gene) & pro_all$class=="enhancer",]
        pro_all$tissue <- gsub('\\..*','',rownames(pro_all))
        #pro_deg <- pro_all[pro_all$gene %in% degs,] #DEG-probe pairs
        if (nrow(pro_all)>0) {
          pro_all$label <- paste0(pro_all$gene,':',pro_all$probe, ':', trait, ':',pro_all$tissue)
          bg <- length(unique(pro_all$label[pro_all$label %in% fdr])) 
          corr <- pro_all[pro_all$p.adj<0.05 & pro_all$label %in% fdr,]
          
          ### gene numbers
          to_plot_2 <- rbind(to_plot_2,c(length(unique(corr[corr$cor>0,"label"]))/bg, "Enhancer", "Positive", trait,length(unique(corr[corr$cor>0,"label"]))),
                             c(length(unique(corr[corr$cor<0,"label"]))/bg, "Enhancer", "Negative", trait,length(unique(corr[corr$cor<0,"label"]))))
          
          bg <- length(unique(pro_all$label))#[!pro_all$label %in% fdr])) 
          corr <- pro_all[pro_all$p.adj<0.05,] # !pro_all$label %in% fdr & 
          ### probe numbers
          to_plot_1 <- rbind(to_plot_1,c(length(unique(corr[corr$cor>0,"label"]))/bg, "Enhancer", "Positive", trait,length(unique(corr[corr$cor>0,"label"]))),
                             c(length(unique(corr[corr$cor<0,"label"]))/bg, "Enhancer", "Negative", trait,length(unique(corr[corr$cor<0,"label"]))))
          
        }
        #Gene Body
        pro_all <- do.call('rbind.data.frame',DMPs_Res[[trait]])
        pro_all <- pro_all[!is.na(pro_all$gene) & pro_all$class=="gene_body",]
        pro_all$tissue <- gsub('\\..*','',rownames(pro_all))
        #pro_deg <- pro_all[pro_all$gene %in% degs,] #DEG-probe pairs
        if (nrow(pro_all)>0) {
          pro_all$label <- paste0(pro_all$gene,':',pro_all$probe, ':', trait, ':',pro_all$tissue)
          bg <- length(unique(pro_all$label[pro_all$label %in% fdr])) 
          corr <- pro_all[pro_all$p.adj<0.05 & pro_all$label %in% fdr,]
          
          ### gene numbers
          to_plot_2 <- rbind(to_plot_2,c(length(unique(corr[corr$cor>0,"label"]))/bg, "Gene Body", "Positive", trait,length(unique(corr[corr$cor>0,"label"]))),
                             c(length(unique(corr[corr$cor<0,"label"]))/bg, "Gene Body", "Negative", trait,length(unique(corr[corr$cor<0,"label"]))))
          
          
          bg <- length(unique(pro_all$label))#[!pro_all$label %in% fdr])) 
          corr <- pro_all[pro_all$p.adj<0.05,] # !pro_all$label %in% fdr & 
          ### probe numbers
          to_plot_1 <- rbind(to_plot_1,c(length(unique(corr[corr$cor>0,"label"]))/bg, "Gene Body", "Positive", trait,length(unique(corr[corr$cor>0,"label"]))),
                             c(length(unique(corr[corr$cor<0,"label"]))/bg, "Gene Body", "Negative", trait,length(unique(corr[corr$cor<0,"label"]))))
          
        }

}



to_plot_1 <- as.data.frame(to_plot_1[c(2:nrow(to_plot_1)),])
to_plot_1$V1 <- 100*as.numeric(to_plot_1$V1)
names(to_plot_1) <- c("N", "type", "Correlation",'Trait',"Number")
to_plot_1$type <- factor(to_plot_1$type, levels = c("Promoter", "Enhancer", "Gene Body"))

#pnominal
library(ggh4x)
traits_cols <- c('#C49122','#4B8C61','#70A0DF','#A76595')
names(traits_cols) <- c("Ancestry", "Sex", "Age", "BMI")
strip <- strip_themed(background_x = elem_list_rect(fill = traits_cols[1:4]))
to_plot_1$Trait <- factor(to_plot_1$Trait, levels = c("Ancestry", "Sex", "Age", "BMI"))
pdf('~/marenostrum/Projects/GTEx_v8/Methylation/Plots/corr_all_n_direction.pnominal.all.pdf', width = 6, height = 2.5)
ggplot(to_plot_1, aes(type, as.numeric(Number), fill=Correlation)) + 
  #geom_col(aes(type, N, fill=Correlation), width = 0.9) +
  geom_bar(stat = 'identity',position = 'fill', alpha=0.8) + 
  xlab("") +
  # ylab("% of DMPs correlated with a DEG") +
  ylab("% of DMPs correlated with a DEG\n in each direction") +
  scale_fill_manual(values=c("#88CCEE", "#CC6677")) + theme_classic() +
  theme(axis.title.y = element_text(margin = margin(r = 2), size = 11),
        axis.text.y = element_text(size = 9, colour = "black"),
        axis.text.x = element_text(size = 9, colour = "black", angle = 90, vjust = 0.5)) +
  geom_text(aes(label=Number, y = as.numeric(Number), x = type), position='fill', stat='identity', size=3,hjust=1, angle=20) +
  #facet_wrap2(~ Trait, strip = strip, nrow = 1, scales = "free_y")
  facet_wrap2(~ Trait, strip = strip, nrow = 1)
dev.off()


to_plot_2 <- as.data.frame(to_plot_2[c(2:nrow(to_plot_2)),])
to_plot_2$V1 <- 100*as.numeric(to_plot_2$V1)
names(to_plot_2) <- c("N", "type", "Correlation",'Trait',"Number")
to_plot_2$type <- factor(to_plot_2$type, levels = c("Promoter", "Enhancer", "Gene Body"))

#FDR
library(ggh4x)
library(ggh4x)
traits_cols <- c('#C49122','#4B8C61','#70A0DF','#A76595')
names(traits_cols) <- c("Ancestry", "Sex", "Age", "BMI")
strip <- strip_themed(background_x = elem_list_rect(fill = traits_cols[1:4]))
to_plot_2$Trait <- factor(to_plot_2$Trait, levels = c("Ancestry", "Sex", "Age", "BMI"))
pdf('~/marenostrum/Projects/GTEx_v8/Methylation/Plots/corr_all_n_direction.fdr.pdf', width = 6, height = 2.5)
ggplot(to_plot_2, aes(type, as.numeric(Number), fill=Correlation)) + 
  #geom_col(aes(type, N, fill=Correlation), width = 0.9) +
  geom_bar(stat = 'identity',position = 'fill', alpha=0.8) + 
  xlab("") +
  # ylab("% of DMPs correlated with a DEG") +
  ylab("% of DMPs correlated with a DEG\n in each direction") +
  scale_fill_manual(values=c("#88CCEE", "#CC6677")) + theme_classic() +
  theme(axis.title.y = element_text(margin = margin(r = 2), size = 11),
        axis.text.y = element_text(size = 9, colour = "black"),
        axis.text.x = element_text(size = 9, colour = "black", angle = 90, vjust = 0.5)) +
  geom_text(aes(label=Number, y = as.numeric(Number), x = type), position='fill', stat='identity', size=3,hjust=1, angle=20) +
  #facet_wrap2(~ Trait, strip = strip, nrow = 1, scales = "free_y")
  facet_wrap2(~ Trait, strip = strip, nrow = 1)
dev.off()



### enrichments ####
library(WebGestaltR)
library(clusterProfiler)
library(org.Hs.eg.db)

genes <- read.delim('~/marenostrum/Projects/GTEx_v8/Methylation/Data/gene_annotation.csv', sep = ',')

results_enrichments <- list()
for (tissue in tissues) {
  for (trait in c("Ancestry", "Sex", "Age", "BMI")) {
    if(tissue %in% sex_tissues & trait == "Sex"){
      next
    } else {
      bg <- unique(DMPs_Res[[trait]][[tissue]]$gene)
      test <- unique(DMPs_Res[[trait]][[tissue]]$gene[DMPs_Res[[trait]][[tissue]]$p.adj<0.05])
      test_down <- enrichGO(gene= test, universe =  bg,
                      keyType = "SYMBOL", OrgDb = org.Hs.eg.db,
                      ont = "BP")
      if (!is.null(test_down)) {
        res <- test_down@result
        results_enrichments[[trait]][[tissue]] <- res[res$p.adjust<0.05,]
      } else{
      results_enrichments[[trait]][[tissue]] <- NA
      }
    }
  }
}

#gsub('\\.*','',genes$ensembl.id[genes$gene.name %in% test])
#gsub('\\.*','',genes$ensembl.id[genes$gene.name %in% bg])

