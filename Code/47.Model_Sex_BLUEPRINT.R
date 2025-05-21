#!/usr/bin/env Rscript
# @Author: Winona Oliveros Diez
# @E-mail: winn95@gmail.com
# @Description: Model sex effects in BLUEPRINT data
# @software version: R=4.2.2

# Parsing
library(optparse)
parser <- OptionParser()
parser <- add_option(parser, opt_str=c("-t", "--tissue"), type="character",
                     dest="tissue",
                     help="Tissue")
options=parse_args(parser)
tissue=options$tissue
# tissue <- "Lung"
tissue <- 'mono'

print(tissue)
first_dir <- "~/2024_GTEx_Methylation/"
#first_dir <- "/gpfs/"
#project_path <- paste0(first_dir, "/projects/bsc83/Projects/GTEx_v8/Methylation/")
project_path <- paste0(first_dir)


print("Reading data")
Sys.time()
data <- read.delim(paste0(project_path,'Data/public/BLUEPRINT/', tissue, "_meth_M_20151028.txt.gz")) #From whole compressed data in 5.6G to compressed 1.4G/1.1Gb only in Lung (the highest number of samples)
Sys.time() #12 minutes to load 15 Gb
M <- data
rownames(M) <- M$meth.id
M$meth.id <- NULL

print("Reading metadata")
metadata <- read.delim(paste0(project_path,'Data/public/BLUEPRINT/', "Blueprint_Epivar_array_data.index"))
metadata <- metadata[metadata$DONOR_ID %in% colnames(M),]
library(dplyr)
metadata <- metadata[,c("DONOR_ID", "DONOR_AGE", "DONOR_SEX")] %>% distinct()

metadata$DONOR_SEX <- as.factor(metadata$DONOR_SEX)
metadata$DONOR_AGE <- as.factor(metadata$DONOR_AGE)

individual_variables <- c("DONOR_SEX", "DONOR_AGE")

rownames(metadata) <- metadata$DONOR_ID
### filter ovary samples ####
metadata$DONOR_ID <- NULL

### make sure order is the same
M <- M[,rownames(metadata)]

metadata_2 <- metadata[,c(individual_variables)]

print("metadata is prepared")

library(limma)
metadata_2$DONOR_AGE <- as.factor(as.numeric(metadata_2$DONOR_AGE))

limma_function <- function(fit, x){
  covariate <<- x #makeContrast does not read the function's environment, so I add covariate to the general environment in my session
  contrast.matrix <- suppressWarnings(makeContrasts(covariate, levels=(fit$design))) #Warnings due to change of name from (Intercept) to Intercept
  fitConstrasts <- suppressWarnings(contrasts.fit(fit, contrast.matrix)) #Warning due to Intercept name
  eb = eBayes(fitConstrasts)
  tt.smart.sv <- topTable(eb,adjust.method = "BH",number=Inf)
  return(tt.smart.sv)
}

mod_2 <- model.matrix( as.formula(paste0("~", paste0(colnames(metadata_2), collapse="+"))), data =  metadata_2)

model_function <- function(mod){
  print("Modelling")
  Sys.time()
  fit_M <- lmFit(M, mod)
  # summary(decideTests(fit_M))
  # to_run <- c("EURv1", "SEX2", "AGE")
  to_run <- c("DONOR_SEXMale")

  print(to_run)
  res_M <- lapply(to_run, function(x) limma_function(fit_M, x) )
  names(res_M) <- to_run
  
  to_return <- res_M
  Sys.time()
  return(to_return)
}

res_2 <- model_function(mod_2) #I would use 5 PEERs

saveRDS(res_2, paste0(project_path, "/Data/Generated_data/", tissue, "_DML_results_Sex_blueprint.rds")) #This is the final model we are using
# saveRDS(res_2, paste0(project_path, "/Tissues/", tissue, "/DML_results_no_BMI.rds")) 
print("Using 5 PEERs:")
print(paste0("  Sex: ", sum(res_2$DONOR_SEXMale$adj.P.Val<0.05)))
print(paste0("  Sex:Hyper: ", nrow(res_2$DONOR_SEXMale[res_2$DONOR_SEXMale$adj.P.Val<0.05 & res_2$DONOR_SEXMale$logFC<0,])))
# print(paste0("  BMI: ", sum(res_2$BMI$adj.P.Val<0.05)))

### enrichment
library(missMethyl)
library(ggplot2)

results <- readRDS('~/2024_GTEx_Methylation/Data/Generated_data/tcel_DML_results_Sex_blueprint.rds')

res <- gometh(rownames(results$DONOR_SEXMale[results$DONOR_SEXMale$adj.P.Val<0.05 & results$DONOR_SEXMale$logFC<0,]), all.cpg=rownames(results$DONOR_SEXMale),
               collection="GO", array.type="EPIC")

res <- res[res$ONTOLOGY=="BP",]
print(table(res$FDR<0.05))
res <- res[res$FDR<0.05,]

# Select top 10 enriched GO terms
top_go <- res %>%
  arrange(FDR) %>%
  slice_head(n = 15) %>%
  mutate(Description = factor(TERM, levels = rev(TERM)))  # to keep order in plot  

## plot
pdf('~/2024_GTEx_Methylation/Data/Generated_data/figure_tcell_up_males_GO_top.pdf', width = 6, height = 4)
ggplot(top_go, aes(x = Description, y = -log10(FDR), fill = DE)) +
  geom_bar(stat = "identity") +
  coord_flip() +
  scale_fill_gradient(low = "lightgreen", high = "darkgreen") +
  labs(title = "Top Enriched GO Terms",
       x = "GO Term",
       y = "-log10(FDR)",
       fill = "Differentially\nMethylated Genes") +
  theme_classic()
dev.off()

### binomial 
binom <- binom.test(2263, 3536, 0.5)
binom <- binom.test(2329, 4017, 0.5)
binom <- binom.test(4864, 13620, 0.5)

