# step3.cell_labels_transferring.R
# --- process scRNAseq data from VMI samples ---
# step 3: Annotate scATAC-seq cells by transferring cell labels from scRNA-seq
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

# load integrated ATAC Seurat object
panc.integrated <- readRDS(file.path(infodir, 'panc.integrated.rds'))

# quantify gene activity
gene.activities <- GeneActivity(panc.integrated, features = VariableFeatures(panc.rna))

# add gene activities as a new assay
panc.integrated[["ACTIVITY"]] <- CreateAssayObject(counts = gene.activities)

# normalize gene activities
DefaultAssay(panc.integrated) <- "ACTIVITY"
panc.integrated %<>% NormalizeData()
panc.integrated %<>% ScaleData(features = rownames(panc.integrated))

# Identify anchors between scRNA-seq and scATAC-seq datasets
transfer.anchors <- FindTransferAnchors(reference = panc.rna, query = panc.integrated, 
                                        features = VariableFeatures(object = panc.rna),
                                        reference.assay = "RNA", query.assay = "ACTIVITY", reduction = "cca")

saveRDS(transfer.anchors, file.path(infodir, 'transfer.anchors.rds'))

# Annotate scATAC-seq cells via label transfer
celltype.predictions <- TransferData(anchorset = transfer.anchors, refdata = panc.rna$merged.cluster,
	weight.reduction = panc.integrated[["integrated_lsi"]], dims = 2:30)

write.table(celltype.predictions, file.path(infodir, 'cell.labels.transferred.txt'), 
            quote=F, sep='\t', col.names=NA)

# add transferred ids and max scores to the seurat object
panc.integrated %<>% AddMetaData(metadata = celltype.predictions %>% 
                                 dplyr::select(c('predicted.id', 'prediction.score.max')))

# save integrated Seurat object
saveRDS(panc.integrated, file.path(infodir, 'panc.integrated.rds'))
