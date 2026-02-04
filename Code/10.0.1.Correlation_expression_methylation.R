#########################################
#### correlation DMG and DEG #####
#########################################
#!/usr/bin/env Rscript
# @Author: Jose Miguel Ramirez & Winona Oliveros
# @E-mail: jose.ramirez1@bsc.es & winona.oliveros@bsc.es
# @Description: Code to correlate gene expression and DNA methylation
# @software version: R=4.2.2

#Loading libraries
library(dplyr)

#Set path 
setwd(system("pwd", intern = T)) #If in linux
# setwd(dirname(rstudioapi::getActiveDocumentContext()$path)) #If using Rstudio
# setwd("..")

#Option 1: Correlate all DMPs to all genes
#Option 2: Correlate all probes to all genes
#Option 3: Correlate all DMPs to all DEGs:
# option <- 3 #To edit depending on the set of probe-genes we want to consider
# print(paste("Option", option))

# Parsing
library(optparse)
parser <- OptionParser()
parser <- add_option(parser, opt_str=c("-t", "--tissue"), type="character",
                     dest="tissue",
                     help="Tissue")
options=parse_args(parser)
tissue=options$tissue
# tissue <- "Lung"

print(tissue)


print("Reading annotation")
#Get a table with probe and a gene
Sys.time()
data_path <- "/gpfs/projects/bsc83/Projects/GTEx_v8/Methylation/Data/"
annotation <- read.delim(paste0(data_path, "Methylation_Epic_gene_promoter_enhancer_processed.txt"), sep = '\t', header = T)
Sys.time()

# Here we have the code to clean the annotation and add nearest gene to unannotated promoters and enhancers
# These edits were saved into the _processed file already
# 
# library(tidyverse)
# library(tidyr)
# library(dplyr)
# annotation <- read.csv("~/marenostrum_scratch/GTEx/v9/Oliva/GPL21145_MethylationEPIC_15073387_v-1-0_processed.csv") ### procesed file
# annotation <- annotation %>% separate_rows(UCSC_RefGene_Name,UCSC_RefGene_Group, sep = ';')
# annotation <- annotation %>% distinct()
# ## I will only consider genes associated to Enhancer/promoter and/or with a gene annotated on the manifest file
# promoter_associated_probes <- annotation[annotation$Regulatory_Feature_Group=='Promoter_Associated',c("IlmnID","UCSC_RefGene_Name","Regulatory_Feature_Group", "Phantom5_Enhancers", "UCSC_RefGene_Group")]
# enhancer_associated_probes <- annotation[annotation$Phantom5_Enhancers!='',c("IlmnID","UCSC_RefGene_Name","Regulatory_Feature_Group", "Phantom5_Enhancers", "UCSC_RefGene_Group")]
# gene_body_associated_probes <- annotation[grep('Body',annotation$UCSC_RefGene_Group),c("IlmnID","UCSC_RefGene_Name","Regulatory_Feature_Group", "Phantom5_Enhancers", "UCSC_RefGene_Group")]
# exon_body_associated_probes <- annotation[grep('1stExon',annotation$UCSC_RefGene_Group),c("IlmnID","UCSC_RefGene_Name","Regulatory_Feature_Group", "Phantom5_Enhancers", "UCSC_RefGene_Group")]
# exon_gene <- rbind(gene_body_associated_probes, exon_body_associated_probes)
# 
# exon_gene <- exon_gene[!exon_gene$IlmnID %in% enhancer_associated_probes$IlmnID,]
# exon_gene <- exon_gene[!exon_gene$IlmnID %in% promoter_associated_probes$IlmnID,]
# enhancer_associated_probes <- enhancer_associated_probes[!enhancer_associated_probes$IlmnID %in% promoter_associated_probes$IlmnID,]
# exon_gene$Type <- 'Gene_Associated'
# enhancer_associated_probes$Type <- 'Enhancer_Associated'
# promoter_associated_probes$Type <- 'Promoter_Associated'
# gene_prom_enh_cpgs <- rbind(exon_gene, enhancer_associated_probes, promoter_associated_probes)
# write.table(gene_prom_enh_cpgs, '~/marenostrum/Projects/GTEx_v8/Methylation/Data/Methylation_Epic_gene_promoter_enhancer_processed.txt', sep = '\t',
#             col.names = T, row.names = F, quote = F)

#From ensembl id to gene symbol
gene_annotation  <- read.delim("/gpfs/projects/bsc83/Projects/GTEx_v8/Laura/00.Data/gencode.v26.GRCh38.genes.biotype_matched_v38.bed")[,c(6,7)]
colnames(gene_annotation) <- c("gene", "symbol")

#These were genes duplicated, I changed their names to their correct one
gene_annotation$symbol[gene_annotation$gene=="ENSG00000253972.5"] <- "MAL2-AS1" #I found this a posteriori
gene_annotation$symbol[gene_annotation$gene=="ENSG00000283992.1"] <- "SLURP2" #Insted of LYNX1
gene_annotation$symbol[gene_annotation$gene=="ENSG00000235271.5"] <- "GC22P027299"
gene_annotation$symbol[gene_annotation$gene=="ENSG00000229694.6"] <- "C9orf73"
gene_annotation$symbol[gene_annotation$gene=="ENSG00000228741.2"] <- "GC13P024553" 


print("Running")
print(tissue)
Sys.time()
results_DML <- readRDS(paste0("Tissues/", tissue, "/DML_results_5_PEERs_continous.rds"))

### read DEG ####
tissues <- c("BreastMammaryTissue", "ColonTransverse" ,"KidneyCortex", "Lung", "MuscleSkeletal" ,"Ovary", "Prostate", "Testis", "WholeBlood")

#GTEx_v8 <- readRDS('/gpfs/projects/bsc83/MN4/bsc83/Projects/ribosomal_proteins/Winona/2022_Ribosomal_analysis/Data/Data_set_1.rds')
GTEx_v8 <- lapply(c(tissues), function(t) readRDS(paste0('/gpfs/projects/bsc83/Projects/GTEx_v8/Methylation/Data/DEA/', t, "/", t, "DEA_results_Ancestry_continous_.results.rds")))
names(GTEx_v8) <- tissues
#head(GTEx_v8)

sexual_tissues <- c("Prostate", "Testis", "Ovary")
if(tissue %in% sexual_tissues){
  print("Sexual tissue")
  individual_variables <- c("EURv1", "AGE", "BMI")
} else{
  individual_variables <- c("EURv1", "SEX2", "AGE", "BMI")
}

traits_meth <- c("EURv1", "SEX2", "AGE", "BMI")
traits_expr <- c("EURv1", "Sex", "Age", "BMI")
names(traits_meth) <- traits_expr
names(traits_expr) <- traits_meth

for (trait in individual_variables) {
  res <- results_DML[[trait]]
  
  signif <- res[res$P.Value<0.05,]
  # table(signif$logFC>0)
  #Reading methylationresiduals
  beta <- readRDS(paste0("/gpfs/projects/bsc83/Projects/GTEx_v8/Methylation/Tissues/", tissue, '/',trait,"_methylation_residuals.continous.rds"))
  Sys.time() 
  
  #Reading expression residuals in the lung:
  expression <- readRDS(paste0("/gpfs/projects/bsc83/Projects/GTEx_v8/Methylation/Tissues//", tissue, '/',traits_expr[trait],"_expression_residuals.continous.rds"))
  
  rownames(expression) <- sapply(rownames(expression), function(gene) gene_annotation$symbol[gene_annotation$gene==gene])
  #From sample id to donor id
  colnames(expression) <- sapply(colnames(expression), function(id) paste0(strsplit(id, "-")[[1]][-3], collapse="-"))
  #Subset expression data to match the donors in DNA methylation data
  expression <- expression[,colnames(expression) %in% colnames(beta)]
  beta <- beta[,colnames(beta) %in% colnames(expression)] 
  identical(colnames(beta), colnames(expression)) #same order
  
  Sys.time() 
  print("Starting to compute correlations")
  
  #Compute correlations
  
  #Now correlate beta values for each position and expression for its associated gene and compute Pearson's correlation
  correlation_function <- function(gene, probe){
    if(!gene %in% rownames(expression)){ return(NA) }
    exp <- expression[rownames(expression)==gene,]
    if(!probe %in% rownames(beta)){ return(NA) }
    bet <- beta[rownames(beta)==probe,]
    test <- cor.test(exp, as.numeric(bet))
    to_return <- list("p.val"= test$p.value, "cor"=test$estimate, "gene"=gene, "probe"=probe)
    return(to_return)
  }
  
  #Running a different analysis per promoter, enhancer and gene body
  annotation_promoter <- annotation[annotation$Type=="Promoter_Associated",]
  annotation_enhancer <- annotation[annotation$Type=="Enhancer_Associated",]
  annotation_gene_body <- annotation[annotation$Type=="Gene_Associated",]
  
  individual_variables <- c("EURv1", "SEX2", "AGE", "BMI")
  if (trait == 'EURv1') {
    trait_deg <- 'Ancestry'
  } 
  if (trait == 'SEX2') {
    trait_deg <- 'Sex'
  } 
  if (trait == 'AGE') {
    trait_deg <- 'Age'
  } 
  if (trait == 'BMI') {
    trait_deg <- 'BMI'
  }
  deg <- GTEx_v8[[tissue]][[trait_deg]][GTEx_v8[[tissue]][[trait_deg]][['adj.P.Val']]<0.05,'gene_name']

  annotation_promoter <- annotation_promoter[annotation_promoter$IlmnID %in% rownames(signif) & annotation_promoter$UCSC_RefGene_Name %in% deg,]
  annotation_enhancer <- annotation_enhancer[annotation_enhancer$IlmnID %in% rownames(signif) & annotation_enhancer$UCSC_RefGene_Name %in% deg,]
  annotation_gene_body <- annotation_gene_body[annotation_gene_body$IlmnID %in% rownames(signif) & annotation_gene_body$UCSC_RefGene_Name %in% deg,]
  #annotation_other <- annotation_other[annotation_other$IlmnID %in% rownames(signif),]

  if (nrow(annotation_promoter)>0) {
    annotation_promoter_df <- annotation_promoter$UCSC_RefGene_Name
    names(annotation_promoter_df) <- annotation_promoter$IlmnID
    annotation_promoter_df <- stack(annotation_promoter_df) #From list to data frame
    annotation_promoter_df <- annotation_promoter_df %>% distinct() #Keep only non duplicated rows
    annotation_promoter_df$ind <- as.character(annotation_promoter_df$ind)
  } else {
    annotation_promoter_df <- NA
  }
  if (nrow(annotation_enhancer)>0) {
    annotation_enhancer_df <- annotation_enhancer$UCSC_RefGene_Name
    names(annotation_enhancer_df) <- annotation_enhancer$IlmnID
    annotation_enhancer_df <- stack(annotation_enhancer_df) #From list to data frame
    annotation_enhancer_df <- annotation_enhancer_df %>% distinct() #Keep only non duplicated rows
    annotation_enhancer_df$ind <- as.character(annotation_enhancer_df$ind)
    
  } else {
    annotation_enhancer_df <- NA
  }
  if (nrow(annotation_gene_body)>0) {
    annotation_gene_body_df <- annotation_gene_body$UCSC_RefGene_Name
    names(annotation_gene_body_df) <- annotation_gene_body$IlmnID
    annotation_gene_body_df <- stack(annotation_gene_body_df) #From list to data frame
    annotation_gene_body_df <- annotation_gene_body_df %>% distinct() #Keep only non duplicated rows
    annotation_gene_body_df$ind <- as.character(annotation_gene_body_df$ind)
  } else {
    annotation_gene_body_df <- NA
  }
  
  if (is.na(annotation_promoter_df)[1]) {
    output_promoters <- list('prom'=list("p.val"= 1, "cor"=0, "gene"=NA, "probe"=NA))
    #names(output_promoters) <- c('prom')
  } else{
    output_promoters <- mapply(function(x,y) correlation_function(x, y), annotation_promoter_df[,1], annotation_promoter_df[,2] , SIMPLIFY = FALSE)
  }
  if (is.na(annotation_enhancer_df)[1]) {
    output_enhancers <- list('enh'=list("p.val"= 1, "cor"=0, "gene"=NA, "probe"=NA))
    #names(output_enhancers) <- c(NA)
  } else{
    output_enhancers <- mapply(function(x,y) correlation_function(x, y), annotation_enhancer_df[,1], annotation_enhancer_df[,2] , SIMPLIFY = FALSE)
  }
  if (is.na(annotation_gene_body_df)[1]) {
    output_gene_body <- list('gene'=list("p.val"= 1, "cor"=0, "gene"=NA, "probe"=NA))
    #names(output_gene_body) <- c(NA)
  } else{
    output_gene_body <- mapply(function(x,y) correlation_function(x, y), annotation_gene_body_df[,1], annotation_gene_body_df[,2] , SIMPLIFY = FALSE)
  }
  #output_other <- mapply(function(x,y) correlation_function(x, y), annotation_other_df[,1], annotation_other_df[,2] )
  
  cleaning_output <- function(output){

    output_2 <- as.data.frame(t(do.call(cbind,output)))
    
    output_2$cor <- unlist(output_2$cor)
    output_2$p.val <- unlist(output_2$p.val)
    output_2$gene <- unlist(output_2$gene)
    output_2$probe <- unlist(output_2$probe)
    
    output_2 <- output_2[!is.na(output_2$p.val),]   #Removing NAs
    
    output_2$p.adj <- p.adjust(output_2$p.val, method = "BH")
    return(output_2)
  }
  output_promoters <- cleaning_output(output_promoters)
  output_promoters$class <- "promoter"
  output_enhancers <- cleaning_output(output_enhancers)
  output_enhancers$class <- "enhancer"
  output_gene_body <- cleaning_output(output_gene_body)
  output_gene_body$class <- "gene_body"
  # output_other <- cleaning_output(output_other)
  # output_other$class <- "other"
  
  print(paste0("Finished computing correlations for", tissue))
  Sys.time()
  
  output <- rbind(output_promoters, output_enhancers, output_gene_body)
  # saveRDS(output, paste0("tissues/Lung/Correlations.rds"))
  
  saveRDS(output, paste0("/gpfs/projects/bsc83/Projects/GTEx_v8/Methylation/Tissues/", tissue, '/',trait,"_Correlations_probes_genes_DEG_DMP.pnominal_ancestry_c.rds"))

}
#---------------------------------------------------










