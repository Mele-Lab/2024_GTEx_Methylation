#!/usr/bin/env Rscript
# @Author: Winona Oliveros Diez
# @E-mail: winn95@gmail.com
# @Description: Plot mediation analysis methylation, demographic traits and gene expression 
# @software version: R=4.2.2

set.seed(1)

# Libraries ----
library(RColorBrewer)
library(ComplexHeatmap)
# GO enrichment --
#library(clusterProfiler)
#library(WebGestaltR)
library(org.Hs.eg.db)
library(ggplot2)
library(cowplot)
# Box plots
library(rcompanion)
library(gtools)

# ---- Data ---- ####

# Individual Traits ----
traits_cols <- c("Age" = "#56B4E9",
                 "Ancestry" = "#E69F00",
                 "Sex" =  "#009E73",
                 "BMI" = "#CC79A7")
traits <- names(traits_cols)

# Tissues ---
first_dir <- "~/cluster/"

project_path <- paste0(first_dir, "Projects/GTEx_v8/Methylation/")

tissues <- c("Lung", "ColonTransverse", "Ovary", "Prostate", "BreastMammaryTissue", "MuscleSkeletal", "KidneyCortex", "Testis", "WholeBlood")
names <- c("Age", "Ancestry", "BMI", "Sex")
#data <- matrix(nrow = length(tissues), ncol=4, dimnames = list(tissues, names))

tissue_info <- readRDS(paste0(first_dir, "Projects/GTEx_v8/Methylation/Data/Tissue_info_whole.rds"))

tissue_info <- tissue_info[!grepl("BreastMammaryTissue_", tissue_info$tissue_ID),]
tissue_info <- tissue_info[tissue_info$tissue_ID %in% tissues,]

# Tissue metadata ----
metadata <- lapply(tissues, function(tissue) readRDS(paste0(project_path,"Tissues/", tissue, "/metadata.rds")))
names(metadata) <- tissues

# ---- Analysis ---- ####
# 1.1 mGene data ----
# inpath_mqtls <- "~/marenostrum_scratch/GTEx/v9/mQTLs/"
# mCpGs <- lapply(tissues, function(tis) {
#   mgene_data <- as.data.frame(data.table::fread(paste0(inpath_mqtls,tis,".mQTLs.conditional.txt.gz")))[mgene_data$V7<0.05 & abs(mgene_data$V3)<250000,]
#   #mgene_data <- mgene_data[mgene_data$V7<0.05 & abs(mgene_data$V3)<250000,]
#   unique(gsub(':.*','',mgene_data$V1[mgene_data$V7<0.05 & abs(mgene_data$V3)<250000]))})
# names(mCpGs) <- tissues


# 1.2 Differential methylation results ----
results_DML <- lapply(tissues, function(tis) 
  readRDS(paste0("~/marenostrum/Projects/GTEx_v8/Methylation/Tissues/",tis,"/DML_results_5_PEERs_continous.rds")))
names(results_DML) <- tissues

# DMP ---
# results_DML_all <- lapply(tissues, function(tis) 
#   results_DML[[tis]][['EURv1']])
# names(results_DML_all) <- tissues
# 
# ## sig results 
# DMP <- lapply(tissues, function(tis) 
#   rownames(results_DML_all[[tis]][results_DML_all[[tis]]$adj.P.Val<0.05,]))
# names(DMP) <- tissues

# 1.3 cis-driven classification ----
inpath <- "~/marenostrum/Projects/GTEx_v8/Methylation/Tissues/"
sex_tissues <- c('Ovary','Prostate','Testis')
for(tissue in tissues){
  #if(!file.exists(paste0(inpath, tissue, "/", tissue, ".Ancestry_DMP.Classification_summary.regmed.rds"))){print(tissue)}
  if(!file.exists(paste0(inpath, tissue, "/", tissue, ".expr_meth.Classification_summary.regmed.pval_ancestry_continous.rds"))){print(tissue)}
}
traits_driven <- c('AGE','EURv1','SEX','BMI')
tissues <- tissues[tissues !='KidneyCortex']
d <- lapply(tissues, function(tissue)
  lapply(traits_driven, function(trait) {
    if (tissue %in% sex_tissues & trait == 'SEX') {
      return(NA)
    } else {
      #readRDS(paste0(inpath, tissue, "/", trait,".expr_meth.Classification_summary.pval005.rds"))
      #readRDS(paste0(inpath, tissue, "/", trait,".expr_meth.Classification_summary.FDR005.rds"))
      readRDS(paste0(inpath, tissue, "/", trait,".expr_meth.Classification_summary.regmed.pval_ancestry_continous.rds"))
    }
  }))

names(d) <- tissues
for (tissue in tissues) {
  names(d[[tissue]]) <- traits
}

# Create dataframe --
data <- unlist(d, recursive = FALSE)
data <- do.call("rbind", data)
data <- as.data.frame(data)
data[,c('Tissue','Trait')] <- stringr::str_split_fixed(rownames(data), '\\.', 2)
data$`Ancestry:DMP` <- as.numeric(data$`Ancestry:DMP`)
data$`Ancestry:DEG` <- as.numeric(data$`Ancestry:DEG`)
data$`Ancestry:DEG:DMP` <- as.numeric(data$`Ancestry:DEG:DMP`)
#data$`Gene:NotModelled` <- as.numeric(data$`Gene:NotModelled`)
data$`Gene:Modelled` <- as.numeric(data$`Gene:Modelled`)

# DMP classified as either cis-driven or not cis-driven --
results <-  lapply(tissues, function(tissue)
  lapply(traits_driven, function(trait) {
    if (tissue %in% sex_tissues & trait == 'SEX') {
      return(NA)
    } else {
      #readRDS(paste0(inpath, tissue, "/", trait,"_DMP.Classified.qval005.rds"))
      #readRDS(paste0(inpath, tissue, "/", trait,"_DMP.Classified.FDR005.rds"))
      readRDS(paste0(inpath, tissue, "/", trait,"_DMP.Classified.regmed.rds"))
    }
  }))

names(results) <- tissues
for (tissue in tissues) {
  names(results[[tissue]]) <- traits
}

res <- unlist(results, recursive = FALSE)
res <- do.call("rbind", res)
res <- as.data.frame(res)
res[,c('Tissue','Trait','gene')] <- stringr::str_split_fixed(rownames(res), '\\.', 3)
res[is.na(res$Class),]

write.table(res, '~/marenostrum/Projects/GTEx_v8/Methylation/Data/mediation_regmed_results.txt', quote = F, col.names = T, row.names = F,
            sep = '\t')

pdf('~/marenostrum/Projects/GTEx_v8/Methylation/Plots/cis_driven_DEGs.FDR.regmed.number_driven.pdf', width = 3, height = 3)
ggplot(res[res$cis_driven>0,], aes(x=cis_driven)) + 
  geom_bar(stat = 'count', colour="black", fill="white")+
  #geom_density(alpha=.2, fill="#FF6666") + scale_x_continuous(breaks=c(1,2,3,4,5,6,7,8,9))
dev.off()

# Expand dataframe --
# data <- cbind.data.frame(data,
#                          sapply(tissues, function(tissue) sapply(traits, function(trait) sum(res$Class=="Cis-driven" & res$Tissue==tissue & res$Trait==trait, na.rm = T))),
#                          sapply(tissues, function(tissue) sapply(traits, function(trait) sum(res$Class=="Not_cis-driven" & res$Tissue==tissue & res$Trait==trait, na.rm = T))),
#                          sapply(tissues, function(tissue) sapply(traits, function(trait) sum(is.na(res$Class)& res$Tissue==tissue & res$Trait==trait, na.rm = T))))
# )
# colnames(data)[10] <- "Cis-driven"
# colnames(data)[11] <- "Not_cis-driven"
# colnames(data)[12] <- "Too_Many_mQTL"
# exprs.data <- data

data[,c("Cis-driven","Not_cis-driven","Too_Many_CpGs")] <- 0
for (tissue in tissues) {
  for (trait in traits) {
    data$`Cis-driven`[data$Tissue==tissue & data$Trait==trait] <- sum(res$Class=="Cis-driven" & res$Tissue==tissue & res$Trait==trait, na.rm = T)
    data$`Not_cis-driven`[data$Tissue==tissue & data$Trait==trait] <- sum(res$Class=="Not_cis-driven" & res$Tissue==tissue & res$Trait==trait, na.rm = T)
    data$Too_Many_CpGs[data$Tissue==tissue & data$Trait==trait] <- sum(is.na(res$Class) & res$Tissue==tissue & res$Trait==trait, na.rm = T)
  }
}

# Ancestry DEG data.frame  ----
dd <- data[,c("Tissue","Trait","Ancestry:DMP",
              "Ancestry:DEG",
              "Ancestry:DEG:DMP",
              "Cis-driven",
              "Not_cis-driven", 
              "Too_Many_CpGs",
              "Gene:Modelled")]
dd[is.na(dd)] <- 0
dd$`Ancestry:DEG:not_DMP` <- as.numeric(dd$`Ancestry:DEG`) - as.numeric(dd$`Ancestry:DEG:DMP`)
dd$perc_cis_driven <- 100*(dd$`Cis-driven`/as.numeric(dd$`Ancestry:DEG:DMP`))
dd$perc_not_cis_driven <- 100*(dd$`Not_cis-driven`/as.numeric(dd$`Ancestry:DEG:DMP`))
round(mean(dd$perc_cis_driven), 2)

DEA_GTEx <- readRDS('~/marenostrum/Projects/GTEx_v8/Methylation/Data/DEA_GTEx.rds')

## number of DEGs
meth_genes <- read.delim('~/marenostrum/Projects/GTEx_v8/Methylation/Data/Methylation_Epic_gene_promoter_enhancer_processed.txt')
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

genes_DE_with_probe_d <- as.data.frame(unlist(genes_DE_with_probe, recursive = T))
genes_DE_with_probe_d[,c('Trait','Tissue')] <- stringr::str_split_fixed(rownames(genes_DE_with_probe_d), '\\.', 2)
colnames(genes_DE_with_probe_d) <- c('#DEGs_with_probe','Trait','Tissue')

dd <- merge(dd, genes_DE_with_probe_d)
dd$perc_cis_driven_2 <- 100*(dd$`Cis-driven`/as.numeric(dd$`#DEGs_with_probe`))

# Data for bar plot --
ddp <- cbind.data.frame("not DEG-DMP" = as.numeric(dd$`Ancestry:DEG`) - as.numeric(dd$`Ancestry:DEG:DMP`), # not sGenes
                        "not cis-driven" = dd$`Not_cis-driven`,
                        "cis-driven" = dd$`Cis-driven`,
                        "not classified" = dd$`Too_Many_CpGs`,
                        "Tissue" = dd$Tissue,
                        "Trait" = dd$Trait# ot modelled: No isQTL with MAF 001, no isQTL with var, no isQTL with no dependance
) 


### dotplot ####
head(dd)
dd$perc_cis_driven_2 <- 100*(dd$`Cis-driven`/(dd$`Cis-driven`+dd$`Not_cis-driven`))
library(reshape2)

tissues_cols <- tissue_info[, 3]
names(tissues_cols) <- tissue_info[, 1]
tissues_cols <- tissues_cols[tissues]
#names(tissues_cols) <- tissue_info$tissue_abbrv

pdf('~/marenostrum/Projects/GTEx_v8/Methylation/Plots/cis_driven_DEGs.FDR.regmed.v2.pdf', width = 4.5, height = 2)
ggplot(data = dd,
       aes(x = Trait,
           y = perc_cis_driven_2,
           color =  Tissue),
) +geom_boxplot(col = "black",
                fill = "white",
                outlier.shape = NA,
                width = 0.8) +
  geom_jitter(alpha = 1,
              size = 2) +
  xlab("") +
  ylab("% cis-driven DEGs") +
  ylim(c(0,40))+
  scale_fill_manual(values = tissues_cols) +
  scale_color_manual(values = tissues_cols) +
  theme(panel.background = element_blank(),
        panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(),
        panel.border = element_rect(colour = "black", fill=NA, size=0.5),
        axis.text = element_text(size = 12),
        axis.title = element_text(size = 14),
        axis.title.y = element_text(size=12),
        plot.title = element_text(hjust = 0.5,
                                  size = 15),
        strip.background = element_rect(fill="#B3B3B3"),
        legend.position = "none") 
dev.off()

pdf('~/marenostrum/Projects/GTEx_v8/Methylation/Plots/cis_driven_DEGs.number.FDR.regmed.pdf', width = 4.5, height = 2)
ggplot(data = dd,
       aes(x = Trait,
           y = `Cis-driven`,
           color =  Tissue),
) +geom_boxplot(col = "black",
                fill = "white",
                outlier.shape = NA,
                width = 0.8) +
  geom_jitter(alpha = 1,
              size = 2) +
  xlab("") +
  ylab("cis-driven DEGs") +
  #ylim(c(0,100.5))+
  scale_fill_manual(values = tissues_cols) +
  scale_color_manual(values = tissues_cols) +
  theme(panel.background = element_blank(),
        panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(),
        panel.border = element_rect(colour = "black", fill=NA, size=0.5),
        axis.text = element_text(size = 12),
        axis.title = element_text(size = 14),
        axis.title.y = ggtext::element_markdown(size=12),
        plot.title = element_text(hjust = 0.5,
                                  size = 15),
        strip.background = element_rect(fill="#B3B3B3"),
        legend.position = "none") 
dev.off()

### number of mediating DMPs ####

pdf("~/marenostrum/Projects/GTEx_v8/Methylation/Plots/Number_DMPS_driving_DEGs_FDR.pdf",
    width = 5, height = 4)
ggplot(res[res$Class=='Cis-driven',], aes(x=cis_driven, fill=Trait)) + 
  geom_bar(stat = "count", width=0.5) + facet_wrap(~Tissue) + ylab('# DEG') + xlab('cis-driving DMPs') +
  scale_x_continuous(breaks = seq(0, 10, by = 1))
dev.off()



