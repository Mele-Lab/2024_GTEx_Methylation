#DEA SUBSTRUCTURES

Sys.setenv(TZ="Europe/Madrid")
# ---------------------- #
start_time <- Sys.time()

# Libraries ####
suppressMessages(library(edgeR))
suppressMessages(library(limma))
suppressMessages(library(parallel))
suppressMessages(library(car))

# Command line arguments ####
args <- commandArgs(trailingOnly=TRUE)

# Functions ####
source("./DEA_and_DSA.R_functions.R")

gene_annotation  <- read.csv("../00.Data/gencode.v39.annotation.bed")
Y_genes <- gene_annotation[gene_annotation$chr=="chrY",]$ensembl.id

# Tissue ----
tissues <- c("BreastMammaryTissue", "ColonTransverse" ,"KidneyCortex", "Lung", "MuscleSkeletal" ,"Ovary", "Prostate", "Testis", "WholeBlood")


for (tissue in tissues){
  outpath <-paste0("/gpfs/projects/bsc83/Projects/GTEx_v8/Methylation/Data/DEA/", tissue, "/")
  print(paste0("---------------------Analyzing ", tissue, "---------------------"))
  if(!dir.exists(outpath)){dir.create(outpath, recursive = T)}
  
  # 1.1 Read in data ----
  #FOR GTEx v8
  counts <- readRDS(paste0("/gpfs/projects/bsc83/Projects/GTEx_v8/Laura/00.Data/v8/1.Final_tissues_to_use/", tissue,"/counts.rds"))
  tpm <- readRDS(paste0("/gpfs/projects/bsc83/Projects/GTEx_v8/Laura/00.Data/v8/1.Final_tissues_to_use/", tissue,"/tpm.rds"))
  metadata <-  readRDS(paste0("/gpfs/projects/bsc83/Projects/GTEx_v8/Laura/00.Data/v8/1.Final_tissues_to_use/", tissue,"/metadata.rds"))
  admixture_ancestry <- read.table('/gpfs/projects/bsc83/Projects/GTEx_v8/Methylation/admixture_inferred_ancestry.txt')
  colnames(admixture_ancestry) <- c('Donor','AFRv1','Ancestry_continous','inferred_ancestry','AFRv2','EURv2')
  metadata <- merge(metadata, admixture_ancestry[,c("Donor","Ancestry_continous")], by='Donor')
  counts <- counts_tot[, metadata$Sample]
  tpm<-tpm_tot[, metadata$Sample]
  gene_annotation  <- read.delim("/gpfs/projects/bsc83/Projects/GTEx_v8/Laura/00.Data/gencode.v26.GRCh38.genes.biotype_matched_v38.bed")
  Y_genes <- gene_annotation[gene_annotation$chr=="chrY",]$ensembl.id
  
  
  
  # 1.2 Genes expressed per tissue ----
  # 1.2.1 TPM>=1 in at least 20% of the tissue samples
  exprs_genes.tpm <- rownames(tpm)[apply(tpm, 1, function(x) sum(x>=1) ) >= 0.2*ncol(tpm)  ]  
  
  # 1.2.2 Count >=10 in at least 20% of the tissue samples
  exprs_genes.counts <- rownames(counts)[ apply(counts, 1, function(x) sum(x>=10) ) >= 0.2*ncol(counts)  ]  
  
  # 1.2.3. Intersect gene lists
  exprs_genes <- intersect(exprs_genes.tpm,
                           exprs_genes.counts) 
  
  # Exclude chrY genes in female-only tissues
  if(tissue %in% c("Uterus","Ovary","Vagina","BreastMammaryTissue_Female")){
    exprs_genes <- exprs_genes[!exprs_genes %in% Y_genes]
  }
  
  # 2 Variables ----
  peers<- c("PEER1", "PEER2")
  tec<- colnames(metadata)[colnames(metadata) %in% c("HardyScale","IschemicTime","RIN","ExonicRate","Cohort", "NucIsoBatch")]
  covariates <- c(peers, tec)
  
  # 2.2 Create DGEList object ----
  dge <- DGEList(counts[exprs_genes,])
  
  # 2.3 Calculate normalization factors (does not do the normalization, only computes the factors) ----
  dge <- calcNormFactors(dge)
  
  # 2.4 Voom ----
  v <- voom(dge, design = NULL, normalize="quantile", save.plot=F, plot = F) # samples are treated as replicates
  
  
  # 4. Differential expression analysis ####
  print("# ---- Running differential expression analysis ---- #")
  
  # 4.1 limma fit : expression ~ covariates + traits ----
  my_data <- list()
  resu<- list()
  summary_results <- list()

  individual_traits <- c("Age","Ancestry_continous", "BMI", "Sex")
  covariates <- c(peers, tec)
  
  fml_args_mod <- paste(c(covariates, individual_traits), collapse = " + ")
  mod <- model.matrix( as.formula(paste(" ~  ", paste(fml_args_mod,collapse = " "))), data =  metadata)
  model<-as.formula(paste(" ~  ", paste0(fml_args_mod,collapse = " ")))
  print(paste0("# ----Model: ",model, " ---- #" ))
  
  # 4.2 Limma fit ----
  fit <- lmFit(v, mod)
  #If we want to save residuals:
  
  # exprs_residuals <- resid(fit, v)
  # saveRDS(exprs_residuals, paste0(outpath,tissue,"_", test, "_residuals.results.rds"))
  
  
  # Add objects to data
  my_data[["dge"]] <- dge
  my_data[["v"]] <- v
  my_data[["fit"]] <- fit
  
  # 4.3 Limma test with interaction term ----
  dea_res <- list()
  
  if (test =="all_cov_inter"){
    
    trait_res <- lapply(c(individual_traits, substructures), function(phenotype) limma_lm(fit, phenotype, metadata))
    names(trait_res) <- c(individual_traits, substructures)
    dea_res <- c(dea_res, trait_res)
    interaction_res <- limma_lm(fit, covariate = NULL, covariate_data = metadata, interaction = paste0(sub, ":Age"))
    name_inter<-paste0("Interaction_Age_", sub)
    dea_res[[name_inter]] <- interaction_res
    new_traits<- c(individual_traits, substructures, name_inter)
    
    
  }else if(test=="one_cov_inter"){
    trait_res <- lapply(c(individual_traits, sub), function(phenotype) limma_lm(fit, phenotype, metadata))
    names(trait_res) <- c(individual_traits, sub)
    dea_res <- c(dea_res, trait_res)
    interaction_res <- limma_lm(fit, covariate = NULL, covariate_data = metadata, interaction = paste0(sub, ":Age"))
    name_inter<-paste0("Interaction_Age_", sub)
    dea_res[[name_inter]] <- interaction_res
    new_traits<- c(individual_traits, sub, name_inter)
    
  }else if(test =="all_cov_no_inter"){
    trait_res <- lapply(c(individual_traits, substructures), function(phenotype) limma_lm(fit, phenotype, metadata))
    names(trait_res) <- c(individual_traits, substructures)
    dea_res <- c(dea_res, trait_res)
    new_traits<- c(individual_traits, substructures)
    
    
  }else if(test =="one_cov_no_inter"){
    trait_res <- lapply(c(individual_traits, sub), function(phenotype) limma_lm(fit, phenotype, metadata))
    names(trait_res) <- c(individual_traits, sub)
    dea_res <- c(dea_res, trait_res)
    new_traits<- c(individual_traits, sub)
    
  }
  
  
  # 5. Computing avrg TPM and exprs var ####
  print("# ---- Calculating avrg TPM and var ---- #")
  
  # 5.1 Compute average TPM expression and variance for each event
  avrg_TPM <- apply(tpm, 1, function(x) mean(log2(x+1)))
  median_TPM <- apply(tpm, 1, function(x) median(log2(x+1)))
  var_TPM <- apply(tpm, 1, function(x) var(x))
  
  
  # Add AvgExprs & ExprsVar and order data.frame
  for(trait in new_traits){
    print(trait)
    if(length(dea_res[[trait]])>1){
      # Add average expression TPM and variance
      dea_res[[trait]]$AvgTPM <- sapply(rownames(dea_res[[trait]]), function(gene) avrg_TPM[gene])
      dea_res[[trait]]$MedianTPM <- sapply(rownames(dea_res[[trait]]), function(gene) median_TPM[gene])
      dea_res[[trait]]$VarTPM <- sapply(rownames(dea_res[[trait]]), function(gene) var_TPM[gene])
      gene_names <- sapply(rownames(dea_res[[trait]]), function(gene) gene_annotation[gene_annotation$ensembl.id==gene, "gene.name.x"])
      names(gene_names) <- NULL
      dea_res[[trait]][["gene.name.x"]] <- gene_names
    }
  }
  
  if (test =="all_cov_no_inter"){
    saveRDS(dea_res,
            paste0(outpath,tissue,"_", test, "_AGE_covariates_and_traits.results.rds"))
    print(paste0("# ---- saved file: ",tissue,"_", test, "_AGE_covariates_and_traits.results.rds", " ---- #") )
  }else{
    saveRDS(dea_res,
            paste0(outpath,tissue,"_", sub, "_", test, "_covariates_and_traits.results.rds"))
    
    print(paste0("# ---- saved file: ",tissue,"_", sub, "_",test, "_covariates_and_traits.results.rds", " ---- #") )
    
  }


}
