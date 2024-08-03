# step2.sample_integration.R
# --- process scRNAseq data from islets samples ---
# step 2: Integrate samples based on LSI
# Author: Tuo Zhang
# Date: 8/1/2024
# 

library(Signac)
library(Seurat)
library(EnsDb.Hsapiens.v86)
library(tidyverse)
library(patchwork)
library(future)
library(magrittr)
library(RColorBrewer)
library(pheatmap)
library(GenomicRanges)

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

# load Seurat objects
seurat.obj.list <- readRDS(file.path(infodir, 'seurat.obj.list.rds'))

# add information to identify dataset of origin
for (k in 1:nrow(sample.info)){
    pid <- sample.info$Name[k]
    seurat.obj.list[[pid]]$dataset <- pid
    seurat.obj.list[[pid]] %<>% RenameCells(add.cell.id = pid)
}

# compute LSI
seurat.obj.list <- lapply(seurat.obj.list, FUN=my.preprocess, min.cutoff=20)

# find integration anchors
integration.anchors <- FindIntegrationAnchors(
    object.list = seurat.obj.list,
    anchor.features = rownames(seurat.obj.list[['BM0']]),
    reduction = "rlsi",
    dims = 2:50
)

saveRDS(integration.anchors, file.path(infodir, 'integration.anchors.rds'))

# integrate LSI embeddings
panc.integrated <- IntegrateEmbeddings(
    anchorset = integration.anchors,
    reductions = panc[["lsi"]],
    new.reduction.name = "integrated_lsi",
    dims.to.integrate = 1:50
)

# create UMAP using the integrated embeddings
panc.integrated %<>% RunUMAP(reduction = "integrated_lsi", dims = 2:50, 
                             n.components=2, seed.use=42, n.neighbors=30, n.epochs=500)

# save integrated Seurat object
saveRDS(panc.integrated, file.path(infodir, 'panc.integrated.rds'))
