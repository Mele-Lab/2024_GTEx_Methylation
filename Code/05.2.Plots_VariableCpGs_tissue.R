
#!/usr/bin/env Rscript
# @Author: Maria Sopena Rios
# @E-mail: maria.sopena@bsc.es
# @Description: Plot enrichemnt of highly variable CpGs across individuals within a tissue 
# @software version: R=4.2.2

library(ggpubr)

first_dir <- "/Users/mariasopenar/cluster/"

tissues <- c("BreastMammaryTissue", "ColonTransverse" ,"KidneyCortex", "Lung", "MuscleSkeletal" ,"Ovary", "Prostate", "Testis", "WholeBlood")

fisher_results <- lapply(tissues, function(tissue) 
  readRDS(paste0(first_dir, "/Projects/GTEx_v8/Methylation/varPart/", tissue, "_highly_variable_CpGs_enrichment_chromHMM.rds")))


#reuse Jose's scripts to plot enrichment of highly variable CpGs per tissue 
read_data_fisher <- function(variables, data, tissue, n){ #Function to prepare data to plot and compute adjusted p value
  
  odds_ratio <- lapply(variables, function(type) data[[type]][['f']]$estimate)
  adj.P.Val <- p.adjust(sapply(variables, function(type) data[[type]][['f']]$p.value), method = "BH")
  CI_down <- lapply(variables, function(type) data[[type]][['f']]$conf.int[1])
  CI_up <- lapply(variables, function(type) data[[type]][['f']]$conf.int[2])
  sample_size <- lapply(variables, function(type) data[[type]][['m']][1,1])
  #sample_size <- n
  
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
  all$tissue <- tissue
  return(all)
}


fisher_results <- do.call(rbind.data.frame, lapply(tissues, function(tissue) {
  x <- readRDS(paste0(first_dir, "/Projects/GTEx_v8/Methylation/varPart/", tissue, "_highly_variable_CpGs_enrichment_chromHMM.rds"))
  n <- readRDS(paste0(first_dir, "/Projects/GTEx_v8/Methylation/varPart/", tissue, "_highly_variable_CpGs.rds") )        
  fish_table <- read_data_fisher(c("Enh","EnhBiv","Het","Quies","ReprPC","TSS","TssBiv","Tx","ZNF/Rpts"), x, tissue, length(n))
  return(fish_table)
}))


colors_traits <- list('AGE'=c('#3D7CD0','#B4D6F6'),
                      'SEX2'=c('#3B734E','#89AA94'),
                      'EURv1'=c('#F0AE21','#F9DE8B'))

fisher_results$sig[fisher_results$sig =="Sig"] <- "FDR < 0.05"
fisher_results$sig[fisher_results$sig =="not Sig"] <- "FDR >= 0.05"
fisher_results$sig <- factor(fisher_results$sig, levels=c("FDR >= 0.05", "FDR < 0.05"))
fisher_results$region <- factor(fisher_results$region, levels=rev(c("Enh","EnhBiv","Het","Quies","ReprPC","TSS","TssBiv","Tx","ZNF/Rpts")))
fisher_results$type[fisher_results$type =="hypo"] <- "Hypomethylation"
fisher_results$type[fisher_results$type =="hyper"] <- "Hypermethylation"

plot_fisher_by_type <- function(type){
  g1 <- ggplot(fisher_results[fisher_results$region==type,], aes(x=log2(oddsRatio), y=tissue, alpha=sig)) +
    geom_errorbar(aes(xmin=log2(CI_down), xmax=log2(CI_up)), width=.3,  color="#1b9e78ff") +
    geom_vline(xintercept = 0) +
    #xlim(0,20) + #Only for Lung to show the 0
    geom_point(size=3, color="#1b9e78ff") + ylab('') + theme_bw() +
    #scale_colour_manual(values=colors_traits[[trait]]) +
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
    ggtitle(paste0(type))+
    scale_y_discrete(breaks=tissues)# + xlim(0, 3)
  
    
    g2 <- ggplot(fisher_results[fisher_results$region==type,]) + geom_col(aes(sample_size, tissue), width = 0.6, fill="#1b9e78ff") +
      theme_classic() + xlab("Number of highly variable CpGs") + ylab("") +
      #scale_fill_manual(values=colors_traits[[trait]]) +
      theme(legend.position = "none",
            axis.text.x = element_text(colour="black", size=13),
            axis.text.y=element_blank(),  #remove y axis labels,
            axis.title.x = element_text(size=16)) +
      scale_x_continuous(n.breaks=3)
    
    p <- ggarrange(g1, g2, labels = c("A", "B"),
                   common.legend = TRUE, legend = "right", widths = c(0.8,0.3))
    pdf(file = paste0("/users/mariasopenar/cluster/Projects/GTEx_v8/Methylation/Plots/chromhmm/enrichment_HighVar_CpG_.pdf"), w = 8, h = 4)
    print(p)
    dev.off()
    return(p)
}


plot_fisher_by_type("Enh")
plot_fisher_by_type("TSS")



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


# plot enrichment of highly variable CpG per tissue

saveRDS(high_var_cpgs, paste0("/gpfs/projects/bsc83/Projects/GTEx_v8/Methylation/varPart/", tissue, "_highly_variable_CpGs.rds"))
saveRDS(names(var_cpg), paste0("/gpfs/projects/bsc83/Projects/GTEx_v8/Methylation/varPart/", tissue, "_all_CpGs.rds"))


library(missMethyl)


GOenrichments <- list()
for (tissue in tissues) {
  print(tissue)
  GOenrichments[[tissue]] <- list()
  high_var_cpgs <- readRDS(paste0(first_dir, "/Projects/GTEx_v8/Methylation/varPart/", tissue, "_highly_variable_CpGs.rds"))
  all_cpgs <- readRDS(paste0(first_dir, "/Projects/GTEx_v8/Methylation/varPart/", tissue, "_all_CpGs.rds"))
  res <- missMethyl::gometh(high_var_cpgs, all.cpg=all_cpgs,
                     collection="GO", array.type="EPIC")
  res <- res[res$ONTOLOGY=="BP",]
  print(table(res$FDR<0.05))
      if (sum(res$FDR<0.05) > 0) {
        GOenrichments[[tissue]] <- res[res$FDR<0.05,]
      }
}


GOenrichments <- readRDS(paste0(first_dir, "/Projects/GTEx_v8/Methylation/varPart/All_tissues_highly_variable_CpGs_Functiona_enrichment.rds"))

for(tissue in names(GOenrichments)){GOenrichments[[tissue]]$tissue <- tissue
GOenrichments[[tissue]]$ID <- rownames(GOenrichments[[tissue]])}

go_df <- do.call(rbind.data.frame, GOenrichments) 

# reduce GO enrichments for plotting 

library(rrvgo);library(dplyr);library(forcats);library(scales)

simMatrix <- calculateSimMatrix(go_df$ID,
                                orgdb="org.Hs.eg.db",
                                ont="BP",
                                method="Rel")

scores <- setNames(-log10(go_df$P.DE), go_df$ID)
go_reduced <- reduceSimMatrix(simMatrix,
                              scores,
                              threshold=0.88,
                              orgdb="org.Hs.eg.db")

go_reduced_all <- merge(go_df, go_reduced, by.x="TERM", by.y="term")

go_reduced_count <- go_reduced_all  %>% group_by(parentTerm, tissue) %>%
  tally() %>%                          # Count occurrences
  mutate(percentage = (n / sum(n)) * 100)
  
go_reduced_count2 <- go_reduced_count %>%
  group_by(parentTerm) %>%
  mutate(shared_tissues = n_distinct(tissue)) %>%
  ungroup() %>%
  mutate(parentTerm = fct_reorder(parentTerm, shared_tissues, .desc = FALSE))

tissues <- c("Lung", "ColonTransverse", "Ovary", "Prostate", "BreastMammaryTissue", "MuscleSkeletal", "KidneyCortex", "Testis", "WholeBlood")
go_reduced_count2$tissue <- factor(go_reduced_count2$tissue , levels = tissues)

ggplot(go_reduced_count2, aes(x = tissue, y = parentTerm, size = n)) +
  geom_point(shape = 21, stroke = 0.5, fill="#1b9e78ff", color="white") +
  scale_size_continuous(name="Number of terms") +
  theme_classic() +
  theme(
    legend.title = element_text(size=12),
    axis.text.x = element_text(colour="black", size=12, angle=90, hjust=1, vjust=0.5),
    axis.text.y = element_text(colour="black", size=12),
    legend.text = element_text(colour="black", size=13),
    axis.title.x = element_text(size=13),
    legend.spacing.y = unit(-0.05, "cm"),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    panel.border = element_rect(colour="black", linewidth=1)
  ) +scale_y_discrete(labels = label_wrap(60))+
  xlab("") + ylab("Parent Terms") +
  labs(title="Highly variable CpGs across individuals")
  


counts <- go_df[go_df$FDR < 0.05,] %>% group_by(TERM) %>%
  tally() %>%                          # Count occurrences
  mutate(percentage = (n / sum(n)) * 100)

counts$Var2 <- ""

ggplot(counts, aes(x = Var2, fill = as.factor(n))) +
  geom_bar() + ylab('Nº of Terms') + xlab('') +
  geom_text(aes(label=after_stat(count), y = after_stat(count)), stat='count', position='stack') +
  labs(fill='Nº of Tissues') + scale_fill_grey(start = 0.9, end = 0) + theme_bw()



