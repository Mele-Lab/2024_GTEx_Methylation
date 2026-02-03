#!/usr/bin/env Rscript
# @Author: Winona Oliveros Diez
# @E-mail: winn95@gmail.com
# @Description: Model demographic traits per tissue
# @software version: R=4.2.2

# Parsing
library(optparse)
parser <- OptionParser()
parser <- add_option(parser, opt_str=c("-t", "--tissue"), type="character",
                     dest="tissue",
                     help="Tissue")
parser <- add_option(parser, opt_str=c("-a", "--ancestry"), type="character",
                      dest="ancestry",
                      help="Ancestry continous / categorical")
parser <- add_option(parser, opt_str=c("-r", "--remove_admixed"), type="logical",
                      dest="remove_admixed",
                      help="Remove")
parser <- add_option(parser, opt_str=c("-s", "--smoking_status"), type="logical",
                      dest="smoking_status",
                      help="smoking information")
options=parse_args(parser)
tissue=options$tissue
ancestry=options$ancestry
remove_admixed <- options$remove_admixed
smoking <- options$smoking_status


# tissue <- "Lung"

print(tissue)
#first_dir <- "~/marenostrum/"
#first_dir <- "/Users/mariasopenar/cluster"
first_dir <- "/gpfs/"
project_path <- paste0(first_dir, "/projects/bsc83/Projects/GTEx_v8/Methylation/")
#project_path <- paste0(first_dir, "/Projects/GTEx_v8/Methylation/")


print("Reading data")
Sys.time()
data <- readRDS(paste0(project_path,'Tissues/', tissue, "/data.rds")) #From whole compressed data in 5.6G to compressed 1.4G/1.1Gb only in Lung (the highest number of samples)
Sys.time() #12 minutes to load 15 Gb
beta <- data

print("Reading metadata")
metadata <- readRDS(paste0(project_path, 'Tissues/', tissue, "/metadata.rds"))

print("Reading Admixture results")
admixture_ancestry <- read.table(paste0(project_path, "/admixture_inferred_ancestry.txt"))
colnames(admixture_ancestry) <- c('SUBJID','AFRv1','Ancestry_continous','inferred_ancestry','AFRv2','EURv2')
metadata <- merge(metadata, admixture_ancestry[,c("SUBJID","Ancestry_continous")], by='SUBJID')

print("Smoking status")
smoking_info <- read.table(paste0(project_path, "Donor_IDs_with_smoking_status.txt"), header = T)
colnames(smoking_info) <- c("SUBJID", "SmokerStatus", "Smoking")
metadata <- merge(metadata, smoking_info, by='SUBJID')

#we can exclude admixed individuals (labeled as AMR in GTEx)
if(remove_admixed==T){
  metadata <- metadata[metadata$Ancestry != "AMR", ]
}

#we can model ancestry as continous or categorical 
if(ancestry=="continous"){
 colnames(metadata) <- gsub("Ancestry_continous", "EURv1", colnames(metadata))
}else{
  colnames(metadata) <- gsub("Ancestry", "EURv1",colnames(metadata))
}

metadata$SEX <- as.factor(metadata$SEX)
if(length(levels(metadata$SEX))==1){
  print("Sexual tissue")
  metadata <- metadata[,-which(names(metadata) == "SEX")]
  individual_variables <- c("EURv1", "AGE", "BMI", "TRISCHD", "DTHHRDY")
  # individual_variables <- c("EURv1", "AGE", "TRISCHD", "DTHHRDY")
} else{
  individual_variables <- c("EURv1", "SEX", "AGE", "BMI", "TRISCHD", "DTHHRDY")
  # individual_variables <- c("EURv1", "SEX", "AGE", "TRISCHD", "DTHHRDY")
}

if(smoking){
  individual_variables <- c(individual_variables, "Smoking")
}

metadata$DTHHRDY <- as.factor(metadata$DTHHRDY)
rownames(metadata) <- metadata$SUBJID
metadata$SUBJID <- NULL

### make sure order is the same
beta <- beta[,rownames(metadata)]
probes <- rownames(beta)
metadata_2 <- metadata[,c("PEER1", "PEER2", "PEER3", "PEER4", "PEER5", individual_variables)]

print("metadata is prepared")

library(limma)

limma_function <- function(fit, x){
  covariate <<- x #makeContrast does not read the function's environment, so I add covariate to the general environment in my session
  contrast.matrix <- suppressWarnings(makeContrasts(covariate, levels=fit$design)) #Warnings due to change of name from (Intercept) to Intercept
  fitConstrasts <- suppressWarnings(contrasts.fit(fit, contrast.matrix)) #Warning due to Intercept name
  eb = eBayes(fitConstrasts)
  tt.smart.sv <- topTable(eb,adjust.method = "BH",number=Inf)
  return(tt.smart.sv)
}

beta <- sapply(beta, as.numeric)
rownames(beta) <- probes
M <- log2(beta/(1-beta)) #-> most papers say M is better for differential although beta should be plotted for interpretation

mod_2 <- model.matrix( as.formula(paste0("~", paste0(colnames(metadata_2), collapse="+"))), data =  metadata_2)

model_function <- function(mod){
  print("Modelling")
  Sys.time()
  fit_M <- lmFit(M, mod)

  to_run <- c("EURv1", "SEX2", "AGE", "BMI")

  if(!"SEX" %in% colnames(metadata)){
    to_run <- to_run[!to_run=="SEX2"]
  }
  print(to_run)
  res_M <- lapply(to_run, function(x) limma_function(fit_M, x) )
  names(res_M) <- to_run

  to_return <- res_M
  Sys.time()
  return(to_return)
}

res_2 <- model_function(mod_2) #I would use 5 PEERs

saveRDS(res_2, paste0(project_path, "/Tissues/", tissue, "/DML_results_5_PEERs_Ancestry_",ancestry, "_remove_admixed_", remove_admixed,"_smoking_", smoking,".peer.rds")) #This is the final model we are using
print(paste0("Using 5 PEERs, ",ancestry, " ancestry and remove admixed individuals ", remove_admixed, ":"))
print(paste0("  EURv1: ", sum(res_2$EURv1$adj.P.Val<0.05)))
print(paste0("  Age: ", sum(res_2$AGE$adj.P.Val<0.05)))

