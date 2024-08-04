# step_1.qc.R
# --- process scRNAseq data from VMI samples ---
# step 1: Load UMI counts data and perform cell QC
# Author: Tuo Zhang
# Date: 8/1/2024
# 

library(Seurat)
library(scater)
library(scran)
library(batchelor)
library(tidyverse)
library(scuttle)
library(future)
library(RColorBrewer)
library(pheatmap)
library(R.utils)
library(magrittr)
library(patchwork)

# folders
workdir <- "."
sourcedir <- file.path(workdir, "data")
figdir <- file.path(workdir, "figure")
infodir <- file.path(workdir, "info")

# project name
project <- "liuliu"

# pattern for defining mitochondrial/ribosomal genes
mito.pattern <- "^MT-"
ribo.pattern <- "^RPL|^RPS"

# random seed
rseed <- 98
set.seed(rseed)

# load functions
setwd(workdir)
source("my_functions.R")

# set parallelization in Seurat
plan("multicore", workers=6)
options(future.globals.maxSize=10*1024^3)

# sample info
sample.info <- data.frame(SeqName=c('non-co-RNA','VMI-UN-RNA','VMI-PRO-RNA'),
                          Group=c('non-coculture','VMI_unstimulated','VMI_pro-inflammatory'),
                          Name=c('BUC','BM0','BM1'))
rownames(sample.info) <- sample.info$Name

# load corrected UMI counts table per patient
raw.counts.list <- list()
for (k in 1:nrow(sample.info)){
    pid <- rownames(sample.info)[k]
    sid <- sample.info$SeqName[k]
    raw.counts.list[[k]] <- my.Read10X(file.path(sourcedir, sid, 'filtered_feature_bc_matrix'), pid)
}
names(raw.counts.list) <- rownames(sample.info)

# merge raw UMI counts tables
raw.counts.all <- my.MergeMatrix.v2(raw.counts.list)

# Initialize the Seurat object with the raw (non-normalized data).
panc.initial <- CreateSeuratObject(counts=raw.counts.all, project=project, assay="RNA", min.cells=0, min.features=0, 
	names.field=1, names.delim="_", meta.data=NULL)

# Calculates the mitochondrial/ribosomal genes per cell
panc.initial[["percent.mito"]] <- PercentageFeatureSet(panc.initial, pattern=mito.pattern)
panc.initial[["percent.ribo"]] <- PercentageFeatureSet(panc.initial, pattern=ribo.pattern)

# Add sample condition
tmeta <- data.frame(row.names=rownames(panc.initial@meta.data))
for (tx in colnames(sample.info)){
  tdic <- as.vector(sample.info[,tx])
  names(tdic) <- rownames(sample.info)
  tmeta[,tx] <- as.vector(tdic[as.vector(panc.initial@meta.data[,"orig.ident"])])
}
panc.initial %<>% AddMetaData(metadata=tmeta)

# perform cell filtering
panc.initial %<>% subset(subset=nFeature_RNA > 300 & nFeature_RNA <= 9000 & 
                         nCount_RNA > 600 & nCount_RNA <= 75000 & percent.mito < 10)

# save Seurat object (initial)
saveRDS(panc.initial, file=file.path(infodir, "panc.initial.rds"))
