#!/usr/bin/env Rscript
# @Author: Winona Oliveros Diez
# @E-mail: winn95@gmail.com
# @Description: Plot comparison model with peers and model with peers + estimates
# @software version: R=4.2.2

## read data and create object
tissues <- c('Lung','ColonTransverse','Prostate','BreastMammaryTissue')

#peer_data <- list()

peer_deconv_data <- lapply(tissues, function(x) readRDS(paste0('~/2024_GTEx_Methylation/Data/Generated_data/',x,'_DML_results_deconv_peer_continous.rds')))
names(peer_deconv_data) <- tissues
peer_data <- lapply(tissues, function(x) readRDS(paste0('~/marenostrum/Projects/GTEx_v8/Methylation/Tissues/',x,'/DML_results_5_PEERs_continous.rds')))
names(peer_data) <- tissues

saveRDS(peer_deconv_data,'~/2024_GTEx_Methylation/Data/Generated_data/results_peers_deconv.rds')
saveRDS(peer_data,'~/2024_GTEx_Methylation/Data/Generated_data/results_peers.rds')

### read objects

obj2 <- readRDS('~/2024_GTEx_Methylation/Data/Generated_data/results_peers_deconv.rds')
obj1 <- readRDS('~/2024_GTEx_Methylation/Data/Generated_data/results_peers.rds')

## extract DMPsfunction
get_dmp_counts <- function(obj, trait) {
  sapply(obj, function(tissue_list) {
    if (!is.null(tissue_list[[trait]])) {
      sum(tissue_list[[trait]]$adj.P.Val<0.05)
    } else {
      NA
    }
  })
}

## get overlap function
get_overlap_counts <- function(obj1, obj2, trait) {
  tissues <- intersect(names(obj1), names(obj2))
  sapply(tissues, function(tissue) {
    dmp1 <- rownames(obj1[[tissue]][[trait]][obj1[[tissue]][[trait]]$adj.P.Val<0.05,])
    dmp2 <- rownames(obj2[[tissue]][[trait]][obj2[[tissue]][[trait]]$adj.P.Val<0.05,])
    if (!is.null(dmp1) && !is.null(dmp2)) {
      length(intersect(dmp1, dmp2))/(min(length(dmp1),length(dmp2)))
    } else {
      NA
    }
  })
}

## create matrix function
build_trait_matrix <- function(trait, obj1, obj2) {
  t1 <- get_dmp_counts(obj1, trait)
  t2 <- get_dmp_counts(obj2, trait)
  ov12 <- get_overlap_counts(obj1, obj2, trait)
  
  all_tissues <- unique(c(names(t1), names(t2), names(ov12)))
  mat <- data.frame(
    Tissue = all_tissues,
    Obj1 = t1[all_tissues],
    Obj2 = t2[all_tissues],
    Overlap_1_2 = ov12[all_tissues]
  )
  rownames(mat) <- mat$Tissue
  mat$Tissue <- NULL
  return(as.matrix(mat))
}

## plot function 
library(pheatmap)

# plot_trait_heatmap <- function(trait, obj1, obj2) {
#   mat <- build_trait_matrix(trait, obj1, obj2)
#   pheatmap(mat,
#            main = paste("Trait:", trait),
#            cluster_rows = TRUE,
#            cluster_cols = FALSE,
#            scale = "none",  # or "none" for raw counts
#            color = colorRampPalette(c("white", "blue"))(100))
# }
# 
# plot_trait_heatmap("AGE", obj1, obj2)
# plot_trait_heatmap("EURv1", obj1, obj2)
# plot_trait_heatmap("SEX2", obj1, obj2)

## complex heatmap
mat <- build_trait_matrix("AGE", obj1, obj2)
colnames(mat) <- c('PEERs','PEERs+Estimates','PropOverlap')

library(ComplexHeatmap)
library(RColorBrewer)
suppressPackageStartupMessages(library(circlize))


my_pretty_num_function <- function(n){
  if(n==""){
    return(n)
  } else{
    prettyNum(n, big.mark = ",")
  }
}

pdf('~/2024_GTEx_Methylation/Data/Generated_data/figure_heatmap_compare_age_peer_estimate.pdf', width = 6, height = 4)
Heatmap(mat,
        heatmap_legend_param = list(legend_height = unit(5, "cm"),
                                    grid_width = unit(1, "cm"),
                                    labels_gp=gpar(fontsize=12),
                                    title_gp=gpar(fontsize=12, fontface=2)),
        col= colorRamp2( c(0,1,max(mat[!is.na(mat)])/4,max(mat[!is.na(mat)])/2,max(mat[!is.na(mat)])),
                         brewer.pal(8, "BuPu")[c(1,2,4,5,7)]),
        # col=colorRamp2( c(0,1,3000),
        #                 c("white","#F5F8FC","#1266b5") ),
        na_col = "white",
        cluster_rows = F,
        cluster_columns = F,
        # name = "DE signal",
        name = "#DMPs",
        row_names_side = "left",
        column_names_side = "top",
        column_names_rot =  60,
        column_names_gp = gpar(fontsize = 12),
        column_names_max_height= unit(9, "cm"),
        row_names_gp = gpar(fontsize = 12),
        # left_annotation = row_ha_left,
        cell_fun = function(j, i, x, y, width, height, fill) {
          grid.text(my_pretty_num_function(mat[i, j]), x, y, gp = gpar(fontsize = 12))}
        
)
dev.off()

### sex
mat <- build_trait_matrix("SEX2", obj1, obj2)
colnames(mat) <- c('PEERs','PEERs+Estimates','PropOverlap')
mat <- mat[c(1:2,4),]

library(ComplexHeatmap)
library(RColorBrewer)
suppressPackageStartupMessages(library(circlize))


my_pretty_num_function <- function(n){
  if(n==""){
    return(n)
  } else{
    prettyNum(n, big.mark = ",")
  }
}

pdf('~/2024_GTEx_Methylation/Data/Generated_data/figure_heatmap_compare_sex_peer_estimate.pdf', width = 6, height = 4)
Heatmap(mat,
        heatmap_legend_param = list(legend_height = unit(5, "cm"),
                                    grid_width = unit(1, "cm"),
                                    labels_gp=gpar(fontsize=12),
                                    title_gp=gpar(fontsize=12, fontface=2)),
        col= colorRamp2( c(0,1,max(mat[!is.na(mat)])/4,max(mat[!is.na(mat)])/2,max(mat[!is.na(mat)])),
                         brewer.pal(8, "BuPu")[c(1,2,4,5,7)]),
        # col=colorRamp2( c(0,1,3000),
        #                 c("white","#F5F8FC","#1266b5") ),
        na_col = "white",
        cluster_rows = F,
        cluster_columns = F,
        # name = "DE signal",
        name = "#DMPs",
        row_names_side = "left",
        column_names_side = "top",
        column_names_rot =  60,
        column_names_gp = gpar(fontsize = 12),
        column_names_max_height= unit(9, "cm"),
        row_names_gp = gpar(fontsize = 12),
        # left_annotation = row_ha_left,
        cell_fun = function(j, i, x, y, width, height, fill) {
          grid.text(my_pretty_num_function(mat[i, j]), x, y, gp = gpar(fontsize = 12))}
        
)
dev.off()

## ancestry 
mat <- build_trait_matrix("EURv1", obj1, obj2)
colnames(mat) <- c('PEERs','PEERs+Estimates','PropOverlap')

library(ComplexHeatmap)
library(RColorBrewer)
suppressPackageStartupMessages(library(circlize))


my_pretty_num_function <- function(n){
  if(n==""){
    return(n)
  } else{
    prettyNum(n, big.mark = ",")
  }
}

pdf('~/2024_GTEx_Methylation/Data/Generated_data/figure_heatmap_compare_ancestry_peer_estimate.pdf', width = 6, height = 4)
Heatmap(mat,
        heatmap_legend_param = list(legend_height = unit(5, "cm"),
                                    grid_width = unit(1, "cm"),
                                    labels_gp=gpar(fontsize=12),
                                    title_gp=gpar(fontsize=12, fontface=2)),
        col= colorRamp2( c(0,1,max(mat[!is.na(mat)])/4,max(mat[!is.na(mat)])/2,max(mat[!is.na(mat)])),
                         brewer.pal(8, "BuPu")[c(1,2,4,5,7)]),
        # col=colorRamp2( c(0,1,3000),
        #                 c("white","#F5F8FC","#1266b5") ),
        na_col = "white",
        cluster_rows = F,
        cluster_columns = F,
        # name = "DE signal",
        name = "#DMPs",
        row_names_side = "left",
        column_names_side = "top",
        column_names_rot =  60,
        column_names_gp = gpar(fontsize = 12),
        column_names_max_height= unit(9, "cm"),
        row_names_gp = gpar(fontsize = 12),
        # left_annotation = row_ha_left,
        cell_fun = function(j, i, x, y, width, height, fill) {
          grid.text(my_pretty_num_function(mat[i, j]), x, y, gp = gpar(fontsize = 12))}
        
)
dev.off()



