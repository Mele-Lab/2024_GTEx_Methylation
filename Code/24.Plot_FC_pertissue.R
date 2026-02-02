#!/usr/bin/env Rscript
# @Author: Winona Oliveros Diez
# @E-mail: winn95@gmail.com
# @Description: Plot FC per tissue and demographic trait for supplementary + plot outliers
# @software version: R=4.2.2

#### plot FC per trait per tissue #####

first_dir <- "marenostrum/"

project_path <- paste0(first_dir, "Projects/GTEx_v8/Methylation/")

project_path <- paste0(basepath, "Projects/GTEx_v8/Methylation/")


tissues <- c("Lung", "ColonTransverse", "Ovary", "Prostate", "BreastMammaryTissue", "MuscleSkeletal", "KidneyCortex", "Testis", "WholeBlood")
names <- c("Age", "Ancestry")
traits_to_use <- c('EURv1','AGE','SEX2')

results_DML <- lapply(tissues, function(tis) 
  readRDS(paste0(project_path,"/Tissues/",tis,"/DML_results_5_PEERs_continous.rds")))
names(results_DML) <- tissues

final_df <- list()

for (tiss in tissues) {
  data <- results_DML[[tiss]]
  data <- do.call(rbind.data.frame, data)
  data$trait <- gsub('\\..*','', rownames(data))
  data <- data[data$adj.P.Val<0.05,]
  final_df[[tiss]] <- data
  
}

final_df <- do.call(rbind.data.frame, final_df)
final_df$tissue <- gsub('\\..*','', rownames(final_df))
head(final_df)

colors_traits <- list('AGE'=c('#3D7CD0'),
                      'SEX2'=c('#3B734E'),
                      'EURv1'=c('#F0AE21'))

### plots 
library(ggridges)
library(ggplot2)
# pdf('plot_test.pdf')
# ggplot(final_df[final_df$trait != 'BMI',], aes(x = abs(logFC), y = trait, fill=trait, height = stat(density))) +
#   geom_density_ridges(rel_min_height = 0.1,stat = "density", scale = 3) + 
#   facet_grid(tissue ~ .) +
#   #xlim(c(-3,3)) + 
#   theme_classic() + 
#   theme(legend.position = 'right')+
#   scale_fill_manual(values = colors_traits)+
#   ylab('')#650x570
# dev.off()

#pdf('marenostrum/Projects/GTEx_v8/Methylation/Plots/Distribution_FC_per_tissue.pdf', width = 12, height = 2.5)
pdf(paste0(basepath, '/Projects/GTEx_v8/Methylation/Plots/Distribution_FC_per_tissue.pdf'), width = 12, height = 2.5)

ggplot(final_df[final_df$trait != 'BMI',], aes(abs(logFC), col=trait, 
                                     group=trait)) + 
  stat_ecdf(geom = "step")+ 
  facet_grid(. ~ tissue) + 
  scale_color_manual(values = c('#3D7CD0','#F0AE21','#3B734E'))+
  labs(title="",
       y = "Empirical Cumulative density", x="abs(logFC)")+
  theme_classic()
dev.off()


# Scale AGE effect sizes to 10 years (so it's comparable to binary traits)
final_df$logFC_scaled <- final_df$logFC
final_df$logFC_scaled[final_df$trait == "AGE"] <- final_df$logFC_scaled[final_df$trait == "AGE"] * 10

pdf(paste0(basepath, '/Projects/GTEx_v8/Methylation/Plots/Distribution_FC_per_tissue_Age10y.pdf'), width = 12, height = 2.5)

ggplot(final_df[final_df$trait != 'BMI',], aes(abs(logFC_scaled), col=trait, 
                                               group=trait)) + 
  stat_ecdf(geom = "step")+ 
  facet_grid(. ~ tissue) + 
  scale_color_manual(values = c('#3D7CD0','#F0AE21','#3B734E'))+
  labs(title="",
       y = "Empirical Cumulative density", x="abs(logFC)")+
  theme_classic()
dev.off()


# Use VarPar estimates to compare traits 
sex_tissues <- c('Ovary','Prostate','Testis')

get_tissue_mean_variation <- function(tissue){
  print(tissue)
  
  hier_res <- readRDS(paste0(project_path, "varPart/",tissue, "_var_part.rds"))
  
  if(tissue %in% sex_tissues){
    tissue_methylation_variation_explained <- hier_res[,traits_prev[-2]]
  }else{
    tissue_methylation_variation_explained <- hier_res[,traits_prev]
  }
  dt <- as.data.table(tissue_methylation_variation_explained, keep.rownames = "feature")
  dt <- melt(dt, id.vars="feature", variable.name="trait", value.name="var_explained")
  dt$tissue <- tissue
  return(dt)  
}

traits_prev <- c('EURv1','SEX','AGE','BMI')
tissue_methylation_variation_explained_all <-rbindlist(lapply(tissues, function(tissue) get_tissue_mean_variation(tissue)))


pdf(paste0(basepath, '/Projects/GTEx_v8/Methylation/Plots/Distribution_VarPar_per_tissue.pdf'), width = 12, height = 2.5)

ggplot(tissue_methylation_variation_explained_all[tissue_methylation_variation_explained_all$trait != 'BMI',], aes(var_explained, col=trait, 
                                               group=trait)) + 
  stat_ecdf(geom = "step")+ 
  facet_grid(. ~ tissue) + 
  scale_color_manual(values = c("AGE"='#3D7CD0',"EURv1"='#F0AE21',"SEX"='#3B734E'))+
  labs(title="",
       y = "Empirical Cumulative density", x="Variance Partition")+
  theme_classic()
dev.off()






##### cpg examples of outlier cpgs #####
annotation <- read.delim('~/marenostrum/Projects/GTEx_v8/Methylation/Data/Methylation_Epic_gene_promoter_enhancer_processed.txt')
View(final_df[
      order( final_df$adj.P.Val, final_df$logFC),
  ])

# Ancestry --
cpgname <- "cg12615916"
cpgname <- 'cg08742194'
cpgname <- 'cg04193820' ### multiple studies etnhicty, autoimmunity, hypomethylated AFR, RNF135 gene
betas <- lapply(tissues, function(tissue) readRDS(paste0(project_path,"Tissues/", tissue, "/data.rds"))[cpgname,])
names(betas) <- tissues

metadata <- lapply(tissues, function(tissue) readRDS(paste0(project_path, "Tissues/",tissue, "/metadata.rds")))
names(metadata) <- tissues

admixture_ancestry <- read.table('~/marenostrum_scratch/GTEx/v8/genotype_data/admixture_inferred_ancestry.txt')
colnames(admixture_ancestry) <- c('SUBJID','AFRv1','EURv1','inferred_ancestry','AFRv2','EURv2')

EA_value <- sapply(tissues, function(tissue)
  median(as.numeric(betas[[tissue]][,colnames(betas[[tissue]]) %in% admixture_ancestry$SUBJID[admixture_ancestry$EURv1>=0.5]]))
)
AA_value <- sapply(tissues, function(tissue)
  median(as.numeric(betas[[tissue]][,colnames(betas[[tissue]]) %in% admixture_ancestry$SUBJID[admixture_ancestry$EURv1<0.5]]))
)

data <- cbind.data.frame("Ancestry" = c(rep("EA",length(EA_value)),
                                        rep("AA",length(AA_value))),
                         "value" = c(EA_value, AA_value),
                         "Tissue" = rep(tissues, 2))
data$Ancestry <- factor(data$Ancestry, levels = c("EA", "AA"), order = T)
data$Tissue <- factor(data$Tissue, levels = tissues,
                      order = T)
tissue_cols <- tissue_info[tissues, "colcodes"]
names(tissue_cols) <- as.character(levels(data$Tissue))
data$cpg <- rep(cpgname, nrow(data))
data$cpg <- factor(data$cpg)

p3 <- ggplot(data = data,
             aes(x = Ancestry,
                 y = (value),
                 col = Tissue),
) +
  geom_violin(col = "black") +
  geom_boxplot(col = "black",
               fill = "white",
               outlier.shape = NA,
               notch = T,
               width = 0.25) +
  geom_point(aes(col = Tissue),
             size = 0.8) +
  geom_line(aes(group=Tissue)) +
  xlab("") +
  ylab("Beta (median)") +
  scale_color_manual(values = tissue_cols) +
  labs(title="") +
  facet_grid(~cpg) +
  theme_bw() +
  theme(panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(),
        axis.text = element_text(size = 12),
        axis.title = element_text(size = 14),
        plot.title = element_text(hjust = 0.5,
                                  size = 15),
        strip.background = element_rect(fill=traits_cols["Ancestry"]),
        legend.position = "none") 

p3


# Age --
cpgname <- "cg16867657"
betas <- lapply(tissues, function(tissue) readRDS(paste0(project_path,"Tissues/", tissue, "/data.rds"))[cpgname,])
names(betas) <- tissues

young_value <- sapply(tissues, function(tissue)
  median(as.numeric(betas[[tissue]][,metadata[[tissue]][as.numeric(metadata[[tissue]]$AGE) < 45, "SUBJID"]]))
)
old_value <- sapply(tissues, function(tissue)
  median(as.numeric(betas[[tissue]][,metadata[[tissue]][as.numeric(metadata[[tissue]]$AGE) >= 45, "SUBJID"]]))
)

data <- cbind.data.frame("Age" = c(rep("[20-45)",length(young_value)),
                                   rep("[45-70]",length(old_value))),
                         "value" = c(young_value, old_value),
                         "Tissue" = rep(tissues, 2))
data$Age <- factor(data$Age, levels = c("[20-45)", "[45-70]"), order = T)
data$Tissue <- factor(data$Tissue, levels = tissues,
                      order = T)
tissue_cols <- tissue_info[tissues, "colcodes"]
names(tissue_cols) <- as.character(levels(data$Tissue))
data$cpg <- rep(cpgname, nrow(data))
data$cpg <- factor(data$cpg)

p4 <- ggplot(data = data,
             aes(x = Age,
                 y = (value),
                 col = Tissue),
) +
  geom_violin(col = "black") +
  geom_boxplot(col = "black",
               fill = "white",
               outlier.shape = NA,
               notch = T,
               width = 0.25) +
  geom_point(aes(col = Tissue),
             size = 0.8) +
  geom_line(aes(group=Tissue)) +
  xlab("") +
  ylab("Beta (median)") +
  scale_color_manual(values = tissue_cols) +
  labs(title="") +
  theme_bw() +
  theme(panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(),
        axis.text = element_text(size = 12),
        axis.title = element_text(size = 14),
        plot.title = element_text(hjust = 0.5,
                                  size = 15),
        strip.background = element_rect(fill=traits_cols["Age"]),
        legend.position = "none") +
  facet_grid(~cpg)

p4
library(ggpubr)
pdf(paste0(plot_path, "Figure_2E.tissue_sharing_examples.pdf"),
    width = 2, height = 8)
ggarrange(p3, p4, ncol = 2)
dev.off()

# Sex --
cpgname <- "cg04858776"

betas <- lapply(tissues, function(tissue) readRDS(paste0(project_path,"Tissues/", tissue, "/data.rds"))[cpgname,])
names(betas) <- tissues

metadata <- lapply(tissues, function(tissue) readRDS(paste0(project_path, "Tissues/",tissue, "/metadata.rds")))
names(metadata) <- tissues

admixture_ancestry <- read.table('~/marenostrum_scratch/GTEx/v8/genotype_data/admixture_inferred_ancestry.txt')
colnames(admixture_ancestry) <- c('SUBJID','AFRv1','EURv1','inferred_ancestry','AFRv2','EURv2')

EA_value <- sapply(tissues, function(tissue)
  median(as.numeric(betas[[tissue]][,colnames(betas[[tissue]]) %in% admixture_ancestry$SUBJID[admixture_ancestry$EURv1>=0.5]]))
)
AA_value <- sapply(tissues, function(tissue)
  median(as.numeric(betas[[tissue]][,colnames(betas[[tissue]]) %in% admixture_ancestry$SUBJID[admixture_ancestry$EURv1<0.5]]))
)

data <- cbind.data.frame("Ancestry" = c(rep("EA",length(EA_value)),
                                        rep("AA",length(AA_value))),
                         "value" = c(EA_value, AA_value),
                         "Tissue" = rep(tissues, 2))
data$Ancestry <- factor(data$Ancestry, levels = c("EA", "AA"), order = T)
data$Tissue <- factor(data$Tissue, levels = tissues,
                      order = T)
tissue_cols <- tissue_info[tissues, "colcodes"]
names(tissue_cols) <- as.character(levels(data$Tissue))
data$cpg <- rep(cpgname, nrow(data))
data$cpg <- factor(data$cpg)

p3 <- ggplot(data = data,
             aes(x = Ancestry,
                 y = (value),
                 col = Tissue),
) +
  geom_violin(col = "black") +
  geom_boxplot(col = "black",
               fill = "white",
               outlier.shape = NA,
               notch = T,
               width = 0.25) +
  geom_point(aes(col = Tissue),
             size = 0.8) +
  geom_line(aes(group=Tissue)) +
  xlab("") +
  ylab("Beta (median)") +
  scale_color_manual(values = tissue_cols) +
  labs(title="") +
  facet_grid(~cpg) +
  theme_bw() +
  theme(panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(),
        axis.text = element_text(size = 12),
        axis.title = element_text(size = 14),
        plot.title = element_text(hjust = 0.5,
                                  size = 15),
        strip.background = element_rect(fill=traits_cols["Ancestry"]),
        legend.position = "none") 

p3


##### plot FC DVPs #######
names <- c("Age", "Ancestry", "Sex")

results_DML <- list()
sex_tissues <- c('Ovary','Prostate','Testis')

for (tissue in tissues) {
  if (tissue %in% sex_tissues) {
    traits <- c("Age", "Ancestry")
  } else {
    traits <- c("Age", "Ancestry", "Sex")
  }
  results_DML[[tissue]] <- lapply(traits, function(trait) (readRDS(paste0("Tissues/", tissue, "/DVP_", trait, ".rds"))))
  names(results_DML[[tissue]]) <- traits
}

final_df <- list()

for (tiss in tissues) {
  data <- results_DML[[tiss]]
  data <- do.call(rbind.data.frame, data)
  data$trait <- gsub('\\..*','', rownames(data))
  data <- data[data$Adj.P.Value<0.05,]
  final_df[[tiss]] <- data
  
}

final_df <- do.call(rbind.data.frame, final_df)
final_df$tissue <- gsub('\\..*','', rownames(final_df))
head(final_df)

colors_traits <- list('Age'=c('#3D7CD0'),
                      'Sex'=c('#3B734E'),
                      'Ancestry'=c('#F0AE21'))

### plots 
library(ggridges)
library(ggplot2)
# pdf('plot_test.pdf')
# ggplot(final_df[final_df$trait != 'BMI',], aes(x = abs(logFC), y = trait, fill=trait, height = stat(density))) +
#   geom_density_ridges(rel_min_height = 0.1,stat = "density", scale = 3) + 
#   facet_grid(tissue ~ .) +
#   #xlim(c(-3,3)) + 
#   theme_classic() + 
#   theme(legend.position = 'right')+
#   scale_fill_manual(values = colors_traits)+
#   ylab('')#650x570
# dev.off()

pdf('~/marenostrum/Projects/GTEx_v8/Methylation/Plots/Distribution_FC_per_tissue.DVPs.pdf', width = 12, height = 2.5)
ggplot(final_df, aes(abs(DiffLevene), col=trait, 
                                               group=trait)) + 
  stat_ecdf(geom = "step")+ 
  facet_grid(. ~ tissue) + 
  scale_color_manual(values = c('#3D7CD0','#F0AE21','#3B734E'))+
  labs(title="",
       y = "Empirical Cumulative density", x="abs(logFC)")+
  theme_classic()
dev.off()


#P


