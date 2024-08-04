# step4.coembed_ATAC_RNA.R
# --- process scRNAseq data from VMI samples ---
# step 4: Coembed scRNA-seq and scATAC-seq cells onto the same plot for visualization purpose
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

# load RNA Seurat object
panc.rna <- readRDS(file.path(rna.infodir, 'panc.rds'))

# load MNN corrected expression data
exp.corrected <- readRDS(file.path(rna.infodir, 'variable.gene.exp.MNN.reconstructed.rds'))

# load transfer anchors
transfer.anchors <- readRDS(file.path(workdir, 'info', "transfer.anchors.rds"))

# load integrated ATAC Seurat object
panc.atac <- readRDS(file.path(infodir, 'panc.integrated.rds'))

# switch assay
DefaultAssay(panc.atac) <- "RNA"

# restrict the imputation to variable genes from scRNA-seq
genes.use <- VariableFeatures(panc.rna)
refdata <- exp.corrected[genes.use, ]

# impute scRNA-seq matrix for ATAC cells
imputation <- TransferData(anchorset=transfer.anchors, refdata=refdata, weight.reduction=panc.atac[["integrated_lsi"]], dims = 2:50)
panc.atac[["RNA"]] <- imputation

# subset to variable genes only
panc.rna.corrected <- panc.rna[genes.use,]

# replace RNA expressions by corrected values
panc.rna.corrected[["RNA"]]@data <- exp.corrected

# merge two Seurat objects
coembed.v4 <- merge(x=panc.rna.corrected, y=panc.atac, add.cell.ids = c("RNA", "ATAC"))

# Run PCA and UMAP on this combined object, to visualize the co-embedding of both
# datasets
coembed.v4 <- ScaleData(coembed.v4, features = genes.use, do.scale = FALSE)
coembed.v4 <- RunPCA(coembed.v4, features = genes.use, verbose = FALSE)
coembed.v4 <- RunUMAP(coembed.v4, dims = 1:30)
coembed.v4 <- RunUMAP(coembed.v4, dims = 1:30, n.components=2, seed.use=42, n.neighbors=30, n.epochs=500)

# save object
saveRDS(coembed.v4, file.path(infodir, 'coembed.v4.rds'))
