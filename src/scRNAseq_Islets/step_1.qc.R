# step_1.qc.R
# --- process scRNAseq data from islets samples ---
# step 1: Load UMI counts data and perform cell QC
# Author: Tuo Zhang
# Date: 8/1/2024
# 

library(scran)
library(Seurat)
library(dplyr)
library(magrittr)
library(MAST)
library(future)
library(batchelor)
library(scater)
library(ggrepel)
library(tibble)
library(tidyr)
library(sceasy)

# folders
workdir <- "."
sourcedir <- file.path(workdir, "data")
figdir <- file.path(workdir, "figure")
infodir <- file.path(workdir, "info")

refdir <- file.path(sourcedir, "reference")

# project name
project <- "liuliu"

# pattern for defining mitochondrial/ribosomal/cov2/cvb4 genes
mito.pattern <- "^MT-"
ribo.pattern <- "^RPL|^RPS"

cov2.pattern <- "^CoV2"
cvb4.pattern <- "polyprotein"

# random seed
rseed <- 98
set.seed(rseed)

# load functions
setwd(workdir)
source("my_functions.R")

# set parallelization in Seurat
plan("multicore", workers=6)
# For certain functions, each worker needs access to certain global variables.
# If these are larger than the default limit, you will see errors.
options(future.globals.maxSize=10*1024^3)

# sample info
sample.info <- data.frame(Name=c("control-9-mock","control-9-SARS-CoV-2","control-10-mock","control-10-SARS-CoV-2","control-10-CVB4","control-11-mock","control-11-SARS-CoV-2","control-12-mock","control-12-CVB4-1","control-12-CVB4-2"), 
                          Condition=c('Mock','COV2','Mock','COV2','CVB4','Mock','COV2','Mock','CVB4','CVB4'), 
                          Donor=c('control-9','control-9','control-10','control-10','control-10','control-11','control-11','control-12','control-12','control-12'),
                          Time=c('24h','24h','24h','24h','24h','24h','24h','24h','24h','48h'))
rownames(sample.info) <- c("D1S1","D1S2","D2S1","D2S2","D2S3","D3S1","D3S2","D4S1","D4S2","D4S3")

# load raw UMI counts table per patient
raw.counts.list <- list()
for (k in 1:nrow(sample.info)){
  pid <- rownames(sample.info)[k]
  sid <- sample.info$Name[k]
  raw.counts.list[[k]] <- my.Read10X(file.path(sourcedir, sid, "filtered_feature_bc_matrix"), pid)
}
names(raw.counts.list) <- rownames(sample.info)

# merge raw UMI counts tables
raw.counts.all <- my.MergeMatrix.v2(raw.counts.list)

# Initialize the Seurat object with the raw (non-normalized data).
panc.initial <- CreateSeuratObject(counts=raw.counts.all, project=project, assay="RNA", min.cells=0, min.features=0, 
                                   names.field=1, names.delim="_", meta.data=NULL)

# Calculates the mitochondrial/ribosomal/viral genes per cell
panc.initial[["percent.mito"]] <- PercentageFeatureSet(panc.initial, pattern=mito.pattern)
panc.initial[["percent.ribo"]] <- PercentageFeatureSet(panc.initial, pattern=ribo.pattern)
panc.initial[["percent.cov2"]] <- PercentageFeatureSet(panc.initial, pattern=cov2.pattern)
panc.initial[["percent.cvb4"]] <- PercentageFeatureSet(panc.initial, pattern=cvb4.pattern)

# Add sample condition
tmeta <- data.frame(row.names=rownames(panc.initial@meta.data))
for (tx in colnames(sample.info)){
  tdic <- as.vector(sample.info[,tx])
  names(tdic) <- rownames(sample.info)
  tmeta[,tx] <- as.vector(tdic[as.vector(panc.initial@meta.data[,"orig.ident"])])
}
panc.initial %<>% AddMetaData(metadata=tmeta)

# perform cell filtering
# nGene > 500, nGene <= 6000, nUMI > 1000, nUMI <= 60000, percent.mito < 15%
panc.initial %<>% subset(subset=nFeature_RNA > 500 & nFeature_RNA <= 6000 & nCount_RNA > 1000 & nCount_RNA <= 60000 & percent.mito < 15)

# separate samples for doublet detection
for (sid in c("D1S1","D1S2","D2S1","D2S2","D2S3","D3S1","D3S2","D4S1","D4S2","D4S3")){
        print(paste('processing', sid))
        # subset seurat object
        tobj <- subset(panc.initial, subset=orig.ident %in% c(sid))
        print(paste(ncol(tobj), 'cells', 'detected.'))

        # convert format
        sceasy::convertFormat(tobj, from="seurat", to="anndata", outFile=file.path(infodir, paste(sid, 'h5ad', sep='.')), assay='RNA', main_layer='counts', transfer_layers=c('counts'))
}

# save Seurat object
saveRDS(panc.initial, file.path(infodir, "panc.initial.rds"))
