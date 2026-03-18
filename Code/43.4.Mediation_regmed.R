#!/usr/bin/env Rscript
# @Author: Winona Oliveros Diez
# @E-mail: winn95@gmail.com
# @Description: Perform mediation analysis methylation, demographic traits and gene expression with regmed
# @software version: R=4.2.2

# Parsing
library(optparse)
library(regmed)
parser <- OptionParser()
parser <- add_option(parser, opt_str=c("-t", "--tissue"), type="character",
                     dest="tissue",
                     help="Tissue")
parser <- add_option(parser, opt_str=c("-d", "--trait"), type="character",
                     dest="trait",
                     help="EURv1", "SEX", "AGE", "BMI")
options=parse_args(parser)
tissue=options$tissue
trait <- options$trait
# tissue <- "Lung"

print(tissue)
#first_dir <- "~/marenostrum/"
first_dir <- "/gpfs/"
project_path <- paste0(first_dir, "/projects/bsc83/Projects/GTEx_v8/Methylation/")
#project_path <- paste0(first_dir, "Projects/GTEx_v8/Methylation/")

start_time <- Sys.time()
print("Reading data")

gene_annotation  <- read.delim("/gpfs/projects/bsc83/Projects/GTEx_v8/Laura/00.Data/gencode.v26.GRCh38.genes.biotype_matched_v38.bed")[,c(9,10)]
#gene_annotation <- read.delim("~/marenostrum/MN4/bsc83/Projects/GTEx_v8/Jose/04_Smoking/github/analysis/data/public/gencode.v26.GRCh38.genes.bed", header=F)[,c(6,7)]

colnames(gene_annotation) <- c("gene", "symbol")

#These were genes duplicated, I changed their names to their correct one
gene_annotation$symbol[gene_annotation$gene=="ENSG00000253972.5"] <- "MAL2-AS1" #I found this a posteriori
gene_annotation$symbol[gene_annotation$gene=="ENSG00000283992.1"] <- "SLURP2" #Insted of LYNX1
gene_annotation$symbol[gene_annotation$gene=="ENSG00000235271.5"] <- "GC22P027299"
gene_annotation$symbol[gene_annotation$gene=="ENSG00000229694.6"] <- "C9orf73"
gene_annotation$symbol[gene_annotation$gene=="ENSG00000228741.2"] <- "GC13P024553" 


print("Reading metadata")
metadata <- readRDS(paste0(project_path, 'Tissues/', tissue, "/metadata.rds"))

print("Reading Admixture results")
admixture_ancestry <- read.table(paste0(project_path,'Data/admixture_inferred_ancestry.txt'))
colnames(admixture_ancestry) <- c('SUBJID','AFRv1','EURv1','inferred_ancestry','AFRv2','EURv2')
metadata <- merge(metadata, admixture_ancestry[,c("SUBJID","EURv1")], by='SUBJID')

sexual_tissues <- c("Prostate", "Testis", "Ovary")

metadata$SEX <- as.factor(metadata$SEX)
if(length(levels(metadata$SEX))==1){
  print("Sexual tissue")
  metadata <- metadata[,-which(names(metadata) == "SEX")]
  individual_variables <- c("EURv1", "AGE", "BMI")
} else{
  individual_variables <- c("EURv1", "SEX", "AGE", "BMI")
}
metadata$DTHHRDY <- as.factor(metadata$DTHHRDY)
rownames(metadata) <- metadata$SUBJID
metadata$SUBJID <- NULL
# if(sum(metadata$Ancestry=="AFR")<5){
#   print("Not enough AFRs")
#   metadata <- metadata[, !colnames(metadata) %in% c("Ancestry")]
#   individual_variables <- individual_variables[!individual_variables %in% "Ancestry"]
# }

### make sure order is the same
#beta <- beta[,rownames(metadata)]

#probes <- rownames(beta)

metadata <- metadata[,c(individual_variables)]

print("metadata expression")
metadata_exp <- readRDS(paste0(project_path,"Tissues/", tissue, "/metadata_expression.rds")) #Samples with both expression and methylation
colnames(admixture_ancestry) <- c('Donor','AFRv1','EURv1','inferred_ancestry','AFRv2','EURv2')
metadata_exp <- merge(metadata_exp, admixture_ancestry[,c("Donor","EURv1")], by='Donor')

metadata_exp$Sex <- as.factor(metadata_exp$Sex)
if(tissue %in% sexual_tissues){
  print("Sexual tissue")
  metadata_exp <- metadata_exp[,-which(names(metadata_exp) == "Sex")]
  individual_variables_exp <- c("EURv1", "Age", "BMI")
} else{
  individual_variables_exp <- c("EURv1", "Sex", "Age", "BMI")
}
metadata_exp$HardyScale <- as.factor(metadata_exp$HardyScale)
rownames(metadata_exp) <- metadata_exp$Donor
metadata_exp$Donor <- NULL

metadata_exp <- metadata_exp[rownames(metadata_exp) %in% rownames(metadata),]
metadata <- metadata[rownames(metadata_exp),]
dim(metadata)
dim(metadata_exp)


# print(head(metadata))
# print(head(metadata_exp))

# Model cpg as:
# meth_res ~ ieVariants_in_eGene + traits
# The genotypes of the sVariants are modelled as a factor, where: 0|0 -> 0, 0|1 -> 1; 1|1 ->2

print("metadata is prepared")
library(tidyverse)
library(dplyr)
library(caret)

###### code for cis-driven expressn
#### Function to fit  linear model per gene, w and w/ independent cis-eQTL of the eGene ####

lm.cis_models <- function(g, meth_residuals, expr_residuals, metadata, gene_variants.list, trait){
  
  print(paste0("----  ", g, "  ----"))
  
  # Get gene residual expression across samples ----
  if(is.vector(expr_residuals)){
    res.g <- as.data.frame(expr_residuals)
  }else{
    res.g <- as.data.frame(t(expr_residuals[g,])) 
  }
  #res.g_long <- data.frame(Individual_ID=colnames(res.g), residuals=res.g[1,], row.names=NULL) #residuals wide to long
  # if (ncol(res.g)>1) {
  #   res.g_long <- reshape2::melt(res.g, variable.name="Individual_ID", value.name = "residuals")
  # } else {
  #   res.g_long <- res.g
  #   res.g_long$Individual_ID <- rownames(res.g_long)
  #   colnames(res.g_long) <- c('residuals','Individual_ID')
  # }
  
  # Merge residuals with metadata -> res.md
  metadata$Individual_ID <- rownames(metadata)
  # res.md <- merge(res.g_long, metadata, by = 'Individual_ID')
  # #res.md.filtered <- res.md[,-which(colnames(res.md)=='Individual_ID')]
  # df <- res.md
  # df <- df %>% mutate_if(is.character,as.factor)
  #print(head(df))
  
  print(paste0('Fitting the model with CpGs'))
  # 1. Select gene cpgs
  # 3. Select traits with variance
  # 4. Select cpgs with variance
  # 5. Run models
  
  # 1. Retrieve the eVariants (ieQTLs in gene) ----
  variants.tissue_g <- gene_variants.list[[g]]
  print(paste0("The gene ", g, " has a total of ",length(variants.tissue_g)," DMPs"))
  
  # ## transform methylation to long format 
  # met_long <- reshape2::melt(meth_residuals)
  # colnames(met_long) <- c('DMP','Individual_ID','meth_residual')
  
  # Retrieve the meth of the DMPs
  vg_g <- t(meth_residuals[variants.tissue_g,])
  if(length(variants.tissue_g)==1){
    vg_g <- as.data.frame((vg_g))
    if (ncol(vg_g) == length(metadata$Individual_ID)) {
      vg_g <- t((vg_g))
      colnames(vg_g) <- variants.tissue_g
    } else {
    colnames(vg_g) <- variants.tissue_g
    }
  }else{
    vg_g <- vg_g 
  }
  #print(head(vg_g))
  
  # 2. Select individual variable
  exposure <- metadata[,trait]
  exposure <- as.numeric(exposure)
  expression <- as.vector(t(res.g))
  lambda.grid <- seq(from = 0.6, to = 0.01, by = -0.05)
  
  # 3. Run models ----
  if(ncol(vg_g)>0){
    # Build models and fit lm ----
    fit.grid <- regmed.grid(exposure, vg_g, expression, lambda.grid, frac.lasso = 0.8)
    #Check if there are linearly dependent sVariants ----
    if(ncol(vg_g)>1){
      # cor matrix ----
      df.snps <- vg_g
      #df.snps %>% mutate_if(is.factor, as.numeric) -> df.snps
      correlationMatrix <- cor(df.snps)
      
      # findCorrelation from {caret} package (filter cutoff=0.9 -> out) ----
      #print('Using findCorrelation from caret package to filter linearly dependent sVariants: cutoff>0.9 -> out ')
      highlyCorrelated <- findCorrelation(correlationMatrix, cutoff=0.9, names = T)
      snps.out <- ifelse(length(highlyCorrelated)==0, NA, paste(highlyCorrelated,collapse = ":"))
      
      if(!is.na(snps.out)){
        print(paste0(snps.out, " is linearly dependent and thus excluded"))
        # Exclude correlated
        snps.f.out <- unlist(strsplit(snps.out, split = ":"))
        snps_in <- colnames(vg_g)[!colnames(vg_g)%in%snps.f.out]    
        if(length(snps_in)==0){
          print(paste0(g, " has 0 DMPs with variance"))
          return(NA)
          break
        }
        # Update models
        vg_g <- vg_g[,snps_in]
        if(is.vector(vg_g)){
          vg_g <- as.data.frame(vg_g)
          colnames(vg_g) <- snps_in
        }else{
          vg_g <- vg_g
        }
        fit.grid <- regmed.grid(exposure, vg_g, expression, lambda.grid, frac.lasso = 0.8)
      }
    }else{
      snps.f.out <- c()
      snps.out <- NA
    } 
    # 6. g report & run analysis ----
    fit.best <- regmed.grid.bestfit(fit.grid)
    summary(fit.best)
    edges.any <- regmed.edges(fit.best, type = "any")
    if (!is.null(edges.any$edges)) {
      pdf(paste0(project_path,'Plots/mediation_',g,'_',trait,'_',tissue,'_edges.pdf'), width = 4, height = 3)
      print(plot.regmed.edges(edges.any))
      dev.off()
    }
    
    # Fit lavaan without penalties
    #expression <- as.vector(t(res.g))
    # dat <- regmed.lavaan.dat(exposure, vg_g, expression)
    # mod.any <- regmed.lavaan.model(edges.any, fit.best)
    # fit.lav.any <- lavaan:::sem(model = mod.any, data = dat)
    # summary.lavaan(fit.lav.any)
    
    # report data frame
    g.df <- data.frame(Tissue = tissue,
                       Trait = trait,
                       symbol = g,
                       indv_total = length(unique(metadata$Individual_ID)),
                       num_DMPs= (ncol(vg_g)),
                       #num_DMPs_lost =  (ncol(vg_g))-length(snps.out[!is.na(snps.out)]),
                       num_DMPs_filtered_out =  ifelse(is.na(snps.out), NA, length(snps.f.out))
    )
    rownames(g.df) <- g
    
    res_summary <- summary(fit.best)
    # 8.output ----
    res <- list("report_summary" = g.df,
                "Mediation_res" = res_summary,
                "DMPs_in" = colnames(vg_g)#,
                # "hier.part" = c(r2.traits, r2.eVariants,
                #                 "n.eVariants" = length(snps_in),
                #                 "EA" = table(res.md[res.md$Individual_ID%in%indv_in,]$EURv1>0.75),
                #                 "AA" = table(res.md[res.md$Individual_ID%in%indv_in,]$EURv1<0.75))
    )
    
    return(res)
    
  }else{
    print(paste0(g, " has 0 snps with variance"))
    return(NA)
    break
  }
  
}



print(trait)
#### DMAnalysis
if (trait == 'SEX') {
  trait_dmp <- 'SEX2'
} else{
  trait_dmp <- trait
}
dea_res <- readRDS(paste0(project_path,'Tissues/', tissue,"/DML_results_5_PEERs_continous.rds"))[[trait_dmp]]
tissues <- c("Lung", "ColonTransverse", "Ovary", "Prostate", "BreastMammaryTissue", "MuscleSkeletal", "KidneyCortex", "Testis", "WholeBlood")
GTEx_v8 <- lapply(c(tissues), function(t) readRDS(paste0(project_path,'/Data/DEA/', t, "/", t, "DEA_results_Ancestry_continous_.results.rds")))
names(GTEx_v8) <- tissues
#GTEx_v8 <- readRDS('~/marenostrum/MN4/bsc83/Projects/ribosomal_proteins/Winona/2022_Ribosomal_analysis/Data/Data_set_1.rds')
for(tis in tissues){names(GTEx_v8[[tis]]) <- gsub("Ancestry_continous", "Ancestry", names(GTEx_v8[[tis]]))}

if (trait == 'EURv1') {
  trait_deg <- 'Ancestry'
} 
if (trait == 'SEX') {
  trait_deg <- 'Sex'
} 
if (trait == 'AGE') {
  trait_deg <- 'Age'
} 
if (trait == 'BMI') {
  trait_deg <- 'BMI'
}

signif <- dea_res[dea_res$P.Value<0.05,]
#signif <- dea_res[dea_res$adj.P.Val<0.05,]
deg <- rownames(GTEx_v8[[tissue]][[trait_deg]][GTEx_v8[[tissue]][[trait_deg]][['adj.P.Val']]<0.05,])
deg_symbol <- (GTEx_v8[[tissue]][[trait_deg]][GTEx_v8[[tissue]][[trait_deg]][['adj.P.Val']]<0.05,"gene.name.x"])

Sys.time()
data_path <- "/gpfs/projects/bsc83/Projects/GTEx_v8/Methylation/Data/"
#data_path <- "~/marenostrum/Projects/GTEx_v8/Methylation/Data/"
annotation <- read.delim(paste0(data_path, "Methylation_Epic_gene_promoter_enhancer_processed.txt"), sep = '\t', header = T)
Sys.time()

gene_probe <- annotation[annotation$IlmnID %in% rownames(signif),]
gene_probe <- gene_probe[gene_probe$UCSC_RefGene_Name %in% deg_symbol,]

# Save summary data ----
d <- list()
d[["Ancestry:DMP"]] <- nrow(signif)
d[["Ancestry:DEG"]] <- length(deg)
d[["Ancestry:DEG:DMP"]] <- length(unique(gene_probe$UCSC_RefGene_Name))

# Print to std.out summary data ----
print(paste0("No. of DMPs: ", nrow(signif)))
print(paste0("No. of DEGs:  ", length(deg)))
print(paste0("No. of DEGs with DMP:  ",  length(unique(gene_probe$UCSC_RefGene_Name))))

# Select differentially methylated cpgs with at least 1 independent mQTL ----
#dea_res <- signif[rownames(signif) %in% gene_probe$IlmnID,] # i-mCpG

expressed_genes <- readRDS(paste0(project_path,'Tissues/',tissue,'/expressed_genes.rds'))
expressed_genes_symbol <- gene_annotation$symbol[gene_annotation$gene %in% expressed_genes]

# Create 'dictionary' (named list)  ----
gene_variants.list <- sapply(unique(gene_probe$UCSC_RefGene_Name), function(x){
  variant <- unique(gene_probe[gene_probe$UCSC_RefGene_Name==x,]$IlmnID)
  return(variant)
}, simplify=F)

print('Finished parsing dictionary')

##### run model #######

library(limma)

genes <- expressed_genes_symbol[expressed_genes_symbol %in% names(gene_variants.list)]
betas <- unique(gene_probe[gene_probe$UCSC_RefGene_Name %in% expressed_genes_symbol,]$IlmnID)
# Recover the expression residuals associated with DEGs ----
meth_residuals <- readRDS(paste0(project_path,'Tissues/', tissue,'/',trait_dmp,"_methylation_residuals.continous.rds"))

# Susbet residuals of genes to be modelled ----
if(length(betas)==1){
  meth_residuals <- as.data.frame(t(as.matrix(meth_residuals[betas,])))
  rownames(meth_residuals) <- betas
}else{
  meth_residuals <- meth_residuals[betas,]    
}

## read expression data 
# Recover the expression residuals associated with DEGs ----
expr_residuals <- readRDS(paste0(project_path,'Tissues/', tissue,'/',trait_dmp,"_expression_residuals.continous.rds"))

# Susbet residuals of genes to be modelled ----

expr_residuals <- as.data.frame((as.matrix(expr_residuals[deg[deg %in% expressed_genes],])))
rownames(expr_residuals) <- deg[deg %in% expressed_genes]

rownames(gene_annotation) <- gene_annotation$gene
genes_symbol <- gene_annotation[rownames(expr_residuals),'symbol']
#we have duplicated symbols therefore we will remove one 
expr_residuals$gene_symbol <- genes_symbol
expr_residuals <- expr_residuals[!duplicated(expr_residuals$gene_symbol),]
expr_residuals <- expr_residuals[!is.na(expr_residuals$gene_symbol),]
rownames(expr_residuals) <- expr_residuals$gene_symbol
expr_residuals <- expr_residuals[,-which(colnames(expr_residuals) == "gene_symbol")]

# keep same colnames for the metadata and the residuals
colnames(expr_residuals) <- sub("-[0-9]+$", "",  colnames(expr_residuals) )

# subset methylation residuals to the common set of donors
meth_residuals <- meth_residuals[, rownames(metadata)]
expr_residuals <- expr_residuals[, rownames(metadata_exp)]


dim(expr_residuals)
dim(meth_residuals)
dim(metadata)
dim(metadata_exp)

#identical(metadata$Sample, colnames(exprs_residuals))
if(!identical(rownames(metadata), colnames(meth_residuals))){
  print("The samples in the metadata and the residual data do not coincide")
  q()
}
if(!identical(rownames(metadata_exp), colnames(expr_residuals))){
  print("The donors in the metadata and the expression data do not coincide")
  q()
}

# Is the Ancestry effect on the mCpG cis-driven or not cis-driven ----
# Is the mCpG cis-driven or not cis-driven ----
print('*****************************************************')
print(paste0(tissue, ' has a total of ', length(genes), ' DEGs with DMPs p.val < 0.05'))
print('Applying glm() function to each of those genes')
print('*****************************************************')
cat('\n')
cat('\n')

# lm per event ----
genes <- rownames(expr_residuals)

library(pbmcapply)

ncores <- as.numeric(Sys.getenv("SLURM_CPUS_PER_TASK"))
if(is.na(ncores)) ncores <- parallel::detectCores() - 1

g.lm_models <- pbmclapply(
  genes,
  function(g) lm.cis_models(
    g, meth_residuals, expr_residuals,metadata, gene_variants.list,trait),mc.cores = ncores,  mc.style = "ETA")
names(g.lm_models) <-  genes
# genes that  cannot modelled ----
#d[["Gene:NotModelled"]] <- sum(is.na(g.lm_models))
# genes modelled ----
g.lm_models <- g.lm_models[!is.na(g.lm_models)] # SNP with no variance or no SNP not correlated
d[["Gene:Modelled"]] <- length(g.lm_models)

# Compare modelA vs modelB ----
deg <- names(g.lm_models)
# Anova
cis_driven <- sapply(deg, function(g) sum((g.lm_models[[g]][["Mediation_res"]]$`alpha*beta`)!=0))

# Parse results data ----
results.df <- cbind.data.frame(deg, cis_driven)
results.df$Class <- ifelse(results.df$cis_driven == 0, "Not_cis-driven","Cis-driven")


# Save ancestry-DEG classified ----
saveRDS(results.df,
        paste0(project_path,'Tissues/',tissue,'/',trait,'_DMP.Classified.regmed.pval.rds'))
saveRDS(g.lm_models,
        paste0(project_path,'Tissues/',tissue,'/',trait,'_DMP.g.lm_models.regmed.pval.rds'))


# Report summary ----
report_df <- do.call(rbind.data.frame,
                     lapply(deg, function(gene)
                       g.lm_models[[gene]][["report_summary"]]))
saveRDS(report_df,
        paste0(project_path,'Tissues/',tissue,'/',trait,'.expr_meth.Report_summary.regmed.pval_ancestry_continous.rds'))

# d ----
saveRDS(d,
        paste0(project_path,'Tissues/',tissue,'/',trait,'.expr_meth.Classification_summary.regmed.pval_ancestry_continous.rds'))



# eVariant not in VCF
# saveRDS(eVariants_lost,
#         paste0(project_path, 'Tissues/',tissue, ".mVariants_not_in_vcf.rds"))

# ---------------------- #
end_time <- Sys.time()
# ---------------------- #

# ---------------------- #
print("")
print("# ---- Elapsed time ---- #")
end_time - start_time
print("# ---------------------- #\n")
# ---------------------- #

# ---------------------- #
print("")
print("# ---- Session info ---- #")
sessionInfo()
print("# ---------------------- #\n")
print("")
# ---------------------- #

