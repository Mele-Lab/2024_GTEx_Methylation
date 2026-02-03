# setwd(dirname(rstudioapi::getActiveDocumentContext()$path)) #If using Rstudio
# setwd("..")
first_dir <- "/gpfs/"
setwd(paste0(first_dir, "/projects/bsc83/Projects/GTEx_v8/Methylation/"))

# -------------- #
print(Sys.time())
#-------------- #

tissues <- c("BreastMammaryTissue", "ColonTransverse" ,"KidneyCortex", "Lung", "MuscleSkeletal" ,"Ovary", "Prostate", "Testis", "WholeBlood")
#tissues <- list.dirs("Tissues/", full.names = F)[-1]
# tissues <- tissues[tissues!="KidneyCortex"]
# tissues <- tissues[33:length(tissues)]

shhh <- suppressPackageStartupMessages
shhh(library(optparse))
option_list = list(
  make_option(c("--model"), action="store", default=NA, type='character'))
opt = parse_args(OptionParser(option_list=option_list))

model <- opt$model

#### reading metadata ####
metadata <- lapply(tissues, function(tissue) readRDS(paste0("Tissues/", tissue, "/metadata.rds")))
names(metadata) <- tissues

#### reading beta values ####
beta <- lapply(tissues, function(tissue) readRDS(paste0("Tissues/", tissue, "/data.rds")))
names(beta) <- tissues

for (tissue in tissues) {
  metadata[[tissue]][,'Tissue'] <- tissue
}

### reading admixture and joining metadata #### 
metadata_df <- do.call("rbind",lapply(tissues, function(tissue) metadata[[tissue]][,c("SUBJID", "Tissue","PEER1", "PEER2", "PEER3", "PEER4", "PEER5")]))
head(metadata_df)
rownames(metadata_df) <- paste0(metadata_df$SUBJID,':',metadata_df$Tissue)

### create label with both subj-id pairs ####
for (tissue in tissues) {
  colnames(beta[[tissue]]) <- paste0(colnames(beta[[tissue]]),':',tissue)
}

### merging all data together ####
beta_df <- do.call("cbind",lapply(tissues, function(tissue) beta[[tissue]]))

### performing linear model #### 

beta_df <- beta_df[,rownames(metadata_df)]
probes <- rownames(beta_df)

# Create M values
beta_df <- sapply(beta_df, as.numeric)
rownames(beta_df) <- probes
M <- log2(beta_df/(1-beta_df)) #-> most papers say M is better for differential although beta should be plotted for interpretation

print('Betas prepared')

### run model ####
library('variancePartition')
library('BiocParallel')
# Specify parallel processing parameters
# this is used implicitly by dream() to run in parallel
param = SnowParam(40, "SOCK", progressbar=TRUE)
register(param)

print('create model')
if(model=="subj_tissue"){
  form <- ~ PEER1 + PEER2 + PEER3 + PEER4 + PEER5 + (1 | SUBJID) + (1 | Tissue)
}else if(model =="subj"){
  form <- ~ PEER1 + PEER2 + PEER3 + PEER4 + PEER5 + (1 | SUBJID)
}else if(model =="tissue_interaction"){
  form <- ~ PEER1 + PEER2 + PEER3 + PEER4 + PEER5 + (1 | Tissue)+(1 | Tissue:SUBJID)
}else if(model =="subj_tissue_interaction"){
  form <- ~ PEER1 + PEER2 + PEER3 + PEER4 + PEER5 + (1 | SUBJID) + (1 | Tissue)+(1 | Tissue:SUBJID)
}
print(form)

### run in chuncks
dfs <- split(as.data.frame(M), (seq(nrow(M))-1) %/% 50000) 
for (i in c(1:length(l))) {
  print('variance Partition') 
  print(paste0('chunck ',i))
  
  vp = fitExtractVarPartModel(as.matrix(dfs[[i]]), form, metadata_df)
  pl <- plotVarPart( sortCols(vp))
  
  saveRDS(vp, paste0(i,'_chunck_var_part_', model, '.rds'))
  
  pdf(file = paste0("Plots/",i,"_chunckvar_part_", model, ".pdf"), w = 6, h = 3.5)
  print(pl)
  dev.off()
}
# print('variance Partition') 
# vp = fitExtractVarPartModel(M, form, metadata_df)
# pl <- plotVarPart( sortCols(vp))
# 
# saveRDS(vp, 'var_part.rds')
# 
# pdf(file = paste0("Plots/var_part.pdf"), w = 6, h = 3.5)
# print(pl)
# dev.off()

