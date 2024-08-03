# step_4.cell_clustering.R
# --- process scRNAseq data from islets samples ---
# step 4: cluster cells and infer cell types
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

# Set a random seed
set.seed(98)

# Cluster the cells
my.dims <- 1:50
panc %<>% FindNeighbors(reduction="mnn", dims=my.dims)
panc %<>% FindClusters(resolution=0.5, verbose=T)

# Run UMAP dimensionality reduction
panc %<>% RunUMAP(dims=my.dims, reduction="mnn", n.components=2, seed.use=42, n.neighbors=30, n.epochs=500)

# set cell identity
panc %<>% SetIdent(value="RNA_snn_res.0.5")

# merge clusters
merge.clust <- c(1,0,0,0,2,3,4,5,0,6,0,0,7,8)
names(merge.clust) <- 0:13

final.clust <- as.vector(merge.clust[as.vector(Idents(panc))])
names(final.clust) <- names(Idents(panc))
final.clust <- factor(final.clust, levels=0:8)

# add final clusters to meta data
panc[["merged.cluster"]] <- final.clust

# set final clusters
Idents(panc) <- "merged.cluster"

# save Seurat object
saveRDS(panc, file.path(infodir, 'panc.rds'))
