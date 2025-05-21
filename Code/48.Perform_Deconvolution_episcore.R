#!/usr/bin/env Rscript
# @Author: Winona Oliveros Diez
# @E-mail: winn95@gmail.com
# @Description: Run Episcore to get cell proportions 
# @software version: R=4.2.2

library(EpiSCORE)

tissue <- 'Prostate'
project_path <- '~/marenostrum/Projects/GTEx_v8/Methylation/'

data(lungSS2mca1)
data(ColonRef)
data(BreastRef)
data(ProstateRef)

### build scRNAseq atlas ## only for Lung
# ncpct.v <- summary(factor(celltypeSS2.idx));
# names(ncpct.v) <- celltypeSS2.v;
# print(ncpct.v);
# 
# ## build expression reference
# expref.o <- ConstExpRef(lungSS2mca1.m,celltypeSS2.idx,celltypeSS2.v,markspecTH=rep(3,4));
# 
# ## build DNAm matrix
# refMscm2.m <- ImputeDNAmRef(expref.o$ref$med,db="SCM2",geneID="SYMBOL");
# refMrmap.m <- ImputeDNAmRef(expref.o$ref$med,db="RMAP",geneID="SYMBOL");
# 
# refMmg.m <- ConstMergedDNAmRef(refMscm2.m,refMrmap.m);
# print(dim(refMmg.m));

### Read DNAm data 
data <- readRDS(paste0(project_path,'Tissues/', tissue, "/data.rds"))
probes <- rownames(data)
beta <- sapply(data, as.numeric)
rownames(beta) <- probes

avSIM.m <- EpiSCORE::constAvBetaTSS(beta,type="850k")

### now perform deconvolution
estF.o <- wRPC(avSIM.m,ref=Prostate_Mref.m,useW=TRUE,wth=0.4,maxit=200);

## saveresults 
saveRDS(estF.o, paste0('~/2024_GTEx_Methylation/Data/Generated_data/',tissue,'_deconvolution_episcore.rds'))

