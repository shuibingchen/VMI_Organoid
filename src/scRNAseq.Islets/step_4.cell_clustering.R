# step_4.cell_clustering.R
# --- process scRNAseq data from islets samples ---
# step 4: cluster cells and infer cell types
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

# folders
workdir <- "."
sourcedir <- file.path(workdir, "data")
figdir <- file.path(workdir, "figure")
infodir <- file.path(workdir, "info")
solo.infodir <- file.path(infodir, 'solo')

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

# load Seurat object
panc <- readRDS(file.path(infodir, "panc.rds"))

# Run non-linear dimensional reduction (UMAP)
panc %<>% RunUMAP(dims=1:50, reduction="mnn", n.components=3, seed.use=42, n.neighbors=30, n.epochs=4000)

# Cluster the cells
panc %<>% FindNeighbors(reduction="mnn", dims=1:50)
panc %<>% FindClusters(resolution=seq(0.05,2,by=0.05), verbose=T)

# set cell identity
panc %<>% SetIdent(value="RNA_snn_res.0.2")

# merge cluster
# C0 + C2                ==>  C0 (Beta cells)
# C1 + C4                ==>  C1 (Alpha cells)
# C3 + C11               ==>  C2 (Acinar cells)
# C5                     ==>  C3 (Ductal cells)
# C6                     ==>  C4 (Mesenchymal cells)
# C7                     ==>  C5 (Delta cells)
# C8 + C12               ==>  C6 (PP cells)
# C9                     ==>  C7 (Immunue cells)
# C10                    ==>  C8 (Endothelial cells)
merge.clust <- c(0,1,0,2,1,3,4,5,6,7,8,2,6)
names(merge.clust) <- 0:12

final.clust <- as.vector(merge.clust[as.vector(Idents(panc))])
names(final.clust) <- names(Idents(panc))
final.clust <- factor(final.clust, levels=0:8)

# add final clusters to meta data
panc[["final.clust"]] <- final.clust

# set final clusters
Idents(panc) <- 'final.clust'

# save seurat object
saveRDS(panc, file=file.path(infodir, "panc.rds"))
