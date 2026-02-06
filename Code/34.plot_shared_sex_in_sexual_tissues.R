#!/usr/bin/env Rscript
# @Author: Winona Oliveros Diez
# @E-mail: winn95@gmail.com
# @Description: plot methylation levels of sex-shared DMPs in sex tissues
# @software version: R=4.2.2

### plot methylation levels shared in sex tissues
first_dir <- "/gpfs/"
setwd(paste0(first_dir, "/projects/bsc83/Projects/GTEx_v8/Methylation/"))# Paths ####
### read sharing <- 
sharing <- readRDS('Data/Sharing_DMP.rds')
sharing_sex <- sharing[sharing$trait=='SEX2' & sharing$number>=2,]
table(sharing_sex$dir)

### read betas 
#setwd(paste0('~/marenostrum/', "Projects/GTEx_v8/Methylation/"))#
tissues <- c('Ovary','Testis','Prostate')
#### reading beta values ####
beta <- lapply(tissues, function(tissue) readRDS(paste0("Tissues/", tissue, "/data.rds")))
names(beta) <- tissues

beta_df <- list()
for (tissue in tissues) {
  colnames(beta[[tissue]]) <- paste0(colnames(beta[[tissue]]),':',tissue)
  beta[[tissue]] <- beta[[tissue]][rownames(beta[[tissue]]) %in% sharing_sex$CG[sharing_sex$dir==1],]
  print(head(beta[[tissue]]))
  probes <- rownames(beta[[tissue]])
  beta[[tissue]] <- sapply(beta[[tissue]], as.numeric)
  rownames(beta[[tissue]]) <- probes
  means <- rowMeans(beta[[tissue]])
  beta_df[[tissue]] <- as.data.frame(means)
  beta_df[[tissue]]$cpg <- rownames(beta_df[[tissue]])
  beta_df[[tissue]]$tis <- tissue
  colnames(beta_df[[tissue]]) <- c('beta','cpg','tissue')
}
to_plot <- do.call("rbind.data.frame",beta_df)
saveRDS(to_plot, 'Data/shared_sex_betas_sexual_tissues.rds')

library(ggplot2)
g <- ggplot(to_plot,aes(tissue, beta)) + geom_violin(scale="width", fill="lightgray") + xlab("") + ylab("Beta values") +
  geom_boxplot(col = "black",
               fill = "white",
               outlier.shape = NA,
               notch = T,
               width = 0.25) +
  geom_jitter(col = "black", 
              alpha = 0.1,
              size = 0.8) +
  #geom_text(data=Summary.data ,aes(x = tissue, y = 1.1, label=n), fontface =2, size = 3) +
  theme_bw() + theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1),
                     panel.grid.major = element_blank(),
                     panel.grid.minor = element_blank(),)
ggsave("Plots/methylation_sex_tissues_shared.pdf", g, device = "pdf", width = 5, height = 4)

## run locally 
to_plot <- readRDS('~/marenostrum/Projects/GTEx_v8/Methylation/Data/shared_sex_betas_sexual_tissues.rds')
## now restric to shared cpgs in polycomb regions
sharing <- readRDS('~/marenostrum/Projects/GTEx_v8/Methylation/Data/Sharing_DMP.rds')
shared_cpgs <- readRDS('~/marenostrum/Projects/GTEx_v8/Methylation/Data/chromHMM_shared_9tissues.rds')

to_plot <- to_plot[to_plot$cpg %in% shared_cpgs$name_ann[shared_cpgs$region_chromhmm %in% c('ReprPC')],]

library(dplyr)
df.summary <- to_plot %>%
  group_by(tissue) %>%
  summarise(
    sd = sd(beta, na.rm = TRUE),
    beta = median(beta)
  )
df.summary

my_comparisons <- list( c("Ovary", "Prostate"), c("Ovary", "Testis"))

g <- ggplot(to_plot,aes(tissue, beta)) + geom_violin(scale="width", fill="lightgray") + xlab("") + ylab("Beta values") +
  geom_boxplot(col = "black",
               fill = "white",
               outlier.shape = NA,
               notch = T,
               width = 0.25) +
  geom_jitter(col = "black", 
              alpha = 0.1,
              size = 0.8) +
  geom_line(aes(group = 1), data = df.summary) +
  stat_compare_means(label = "p.format", paired = TRUE, comparisons = my_comparisons) +
  theme_bw() + theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1),
                     panel.grid.major = element_blank(),
                     panel.grid.minor = element_blank(),)
ggsave("~/marenostrum/Projects/GTEx_v8/Methylation/Plots/methylation_sex_tissues_shared.v2.pvalues.pdf", g, device = "pdf", width = 5, height = 4)

library(ggpubr)
ggpaired(to_plot[to_plot$tissue %in% c('Ovary','Testis'),], x = "tissue", y = "beta",
         line.color = "gray", line.size = 0.4, color = 'black', fill='white',
         palette = "jco")+
  stat_compare_means(paired = TRUE)

ggpaired(to_plot[to_plot$tissue %in% c('Ovary','Prostate'),], x = "tissue", y = "beta",
         line.color = "gray", line.size = 0.4, color = 'black', fill='white',
         palette = "jco")+
  stat_compare_means(paired = TRUE)

