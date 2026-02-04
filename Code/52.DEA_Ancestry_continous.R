
#!/usr/bin/env Rscript
# @Author: Maria Sopena Rios
# @E-mail: maria.sopena@bsc.es
# @Description: Run Differential expression analysis with GTEx samples 
# @software version: R=4.2.2


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
#source("./DEA_and_DSA.R_functions.R")
limma_lm <- function(fit, covariate, covariate_data, interaction = NULL) {
  v.contrast <- rep(0, ncol(fit$design))
  
  if (!is.null(interaction)) {
    # If testing for interaction, extract the interaction term
    interaction_term <- interaction
    interaction_idx <- grep(interaction_term, colnames(fit$design))
    if (length(interaction_idx) > 0) {
      v.contrast[interaction_idx] <- 1
      contrast.matrix <- cbind("C1" = v.contrast)
    } else {
      stop("Interaction term not found in the design matrix.")
    }
  } else {
    # Original covariate testing
    if (is.factor(covariate_data[, covariate])) {
      if (covariate == "Sex") {
        if (table(covariate_data$Sex)["1"] >= 10 & table(covariate_data$Sex)["2"] >= 10) {
          v.contrast[which(colnames(fit$design) == "Sex2")] <- 1
          contrast.matrix <- cbind("C1" = v.contrast)
        } else {
          return(NA)
        }
      }else if (covariate == "Age_bin") {
        if (table(covariate_data$Age_bin)["1"] >= 5 & table(covariate_data$Age_bin)["0"] >= 5) {
          v.contrast[which(colnames(fit$design) == "Age_bin1")] <- 1
          contrast.matrix <- cbind("C1" = v.contrast)
        } else {
          return(NA)
        }
        
      }else if (covariate == "Ancestry") {
        if (table(covariate_data$Ancestry)["AFR"] >= 10 & table(covariate_data$Ancestry)["EUR"] >= 10) {
          v.contrast[which(colnames(fit$design) == "AncestryAFR")] <- 1
          contrast.matrix <- cbind("C1" = v.contrast)
        } else {
          return(NA)
        }
      }else if (covariate == "Adenomyosis") {
        if (table(covariate_data$Adenomyosis)["1"] >= 5 & table(covariate_data$Adenomyosis)["0"] >= 5) {
          v.contrast[which(colnames(fit$design) == "Adenomyosis1")] <- 1
          contrast.matrix <- cbind("C1" = v.contrast)
        } else {
          return(NA)
        }
      }else {
        return(NULL)
      }
    } else {
      # Continuous variable
      v.contrast[which(colnames(fit$design) == covariate)] <- 1
      contrast.matrix <- cbind("C1" = v.contrast)
    }
  }
  
  fitContrasts <- contrasts.fit(fit, contrast.matrix)
  eb <- eBayes(fitContrasts)
  tt.smart.sv <- topTable(eb, adjust.method = "BH", number = Inf)
  return(tt.smart.sv)
}


gene_annotation  <- read.delim("/gpfs/projects/bsc83/Projects/GTEx_v8/Laura/00.Data/gencode.v26.GRCh38.genes.biotype_matched_v38.bed")
Y_genes <- gene_annotation[gene_annotation$chr=="chrY",]$ensembl.id

# Tissue ----
tissues <- c("BreastMammaryTissue", "ColonTransverse" ,"KidneyCortex", "Lung", "MuscleSkeletal" ,"Ovary", "Prostate", "Testis", "WholeBlood")
tissues <- c("Testis", "WholeBlood")


for (tissue in tissues){
  outpath <-paste0("/gpfs/projects/bsc83/Projects/GTEx_v8/Methylation/Data/DEA/", tissue, "/")
  print(paste0("---------------------Analyzing ", tissue, "---------------------"))
  if(!dir.exists(outpath)){dir.create(outpath, recursive = T)}
  
  # 1.1 Read in data ----
  #FOR GTEx v8
  print("Reading input data")
  counts <- readRDS(paste0("/gpfs/projects/bsc83/Projects/GTEx_v8/Laura/00.Data/v8/1.Final_tissues_to_use/", tissue,"/counts.rds"))
  tpm <- readRDS(paste0("/gpfs/projects/bsc83/Projects/GTEx_v8/Laura/00.Data/v8/1.Final_tissues_to_use/", tissue,"/tpm.rds"))
  metadata <-  readRDS(paste0("/gpfs/projects/bsc83/Projects/GTEx_v8/Laura/00.Data/v8/1.Final_tissues_to_use/", tissue,"/metadata.rds"))
  admixture_ancestry <- read.table('/gpfs/projects/bsc83/Projects/GTEx_v8/Methylation/admixture_inferred_ancestry.txt')
  colnames(admixture_ancestry) <- c('Donor','AFRv1','Ancestry_continous','inferred_ancestry','AFRv2','EURv2')
  metadata <- merge(metadata, admixture_ancestry[,c("Donor","Ancestry_continous")], by='Donor')
  counts <- counts[, metadata$Sample]
  tpm<-tpm[, metadata$Sample]
  
  
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

  if(! tissue %in% c("Vagina","Uterus","Ovary","Prostate","Testis","BreastMammaryTissue_Female","BreastMammaryTissue_Male")){
    individual_traits <- c("Age","Ancestry_continous", "Sex","BMI")
  }else{
    individual_traits <- c("Age","Ancestry_continous", "BMI")
  }
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
  trait_res <- lapply(c(individual_traits), function(phenotype) limma_lm(fit, phenotype, metadata))
  names(trait_res) <- c(individual_traits)
  dea_res <- c(dea_res, trait_res)
  new_traits<- c(individual_traits)

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
  
saveRDS(dea_res,
        paste0(outpath,tissue,"DEA_results_Ancestry_continous_.results.rds"))
    
    print(paste0("# ---- saved file: ",tissue,"DEA_results_Ancestry_continous_.results.rds", " ---- #") )
    
  }



