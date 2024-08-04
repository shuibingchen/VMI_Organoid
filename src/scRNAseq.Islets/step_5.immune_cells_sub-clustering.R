# step_5.immune_cells_sub-clustering.R
# --- process scRNAseq data from islets samples ---
# step 5: sub-clustering on the immune cell population
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
library(RColorBrewer)

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

# load Seurat object
panc <- readRDS(file.path(infodir, "panc.rds"))

# fetch immune cells
immune.cells.info <- FetchData(panc, vars=c('Name','ident')) %>% rownames_to_column('cellID') %>% filter(ident == 7)

# load rescaled SingleCellExperiment object
rescaled.sce.list <- readRDS(file=file.path(infodir, "rescaled.sce.list.rds"))

# subset immune cells
rescaled.sce.list <- lapply(rescaled.sce.list, FUN=function(x) { x[, intersect(colnames(x), immune.cells.info$cellID)] })

# create a seurat object with raw UMI counts
panc.immune <- CreateSeuratObject(counts=as(do.call(cbind, lapply(rescaled.sce.list, function(x) counts(x))), "dgCMatrix"), 
	project=project, assay="RNA", min.cells=0, min.features=0,
	names.field=1, names.delim="_", meta.data=NULL)

# Calculates the mitochondrial/ribosomal/virus genes per cell
panc.immune[["percent.mito"]] <- PercentageFeatureSet(panc.immune, pattern=mito.pattern)
panc.immune[["percent.ribo"]] <- PercentageFeatureSet(panc.immune, pattern=ribo.pattern)
panc.immune[["percent.cov2"]] <- PercentageFeatureSet(panc.immune, pattern=cov2.pattern)
panc.immune[["percent.cvb4"]] <- PercentageFeatureSet(panc.immune, pattern=cvb4.pattern)

# Add sample condition
tmeta <- data.frame(row.names=rownames(panc.immune@meta.data))
for (tx in colnames(sample.info)){
  tdic <- as.vector(sample.info[,tx])
  names(tdic) <- rownames(sample.info)
  tmeta[,tx] <- as.vector(tdic[as.vector(panc.immune@meta.data[,"orig.ident"])])
}
panc.immune %<>% AddMetaData(metadata=tmeta)

# replace normalized data with the scran normalized data
panc.immune[["RNA"]]@data <- as(do.call(cbind, lapply(rescaled.sce.list, function(x) logcounts(x))) * log(2), "dgCMatrix")

# Identification of highly variable features (feature selection)
panc.immune %<>% FindVariableFeatures(selection.method="vst", nfeatures=3500)

# remove dissociation-related genes and ribosomal genes from variable gene list
# and select the remaining top 3000 genes for MNN-based correction
variable.genes <- setdiff(VariableFeatures(panc.immune), c(disso.genes, grep(mito.pattern, rownames(panc.immune), value=T), 
                                                           grep(ribo.pattern, rownames(panc.immune), value=T), 
                                                           grep(cov2.pattern, rownames(panc.immune), value=T),
                                                           grep(cvb4.pattern, rownames(panc.immune), value=T)))
variable.genes <- head(variable.genes, 3000)

# perform MMN-based correction
original <- lapply(rescaled.sce.list, function(x) {logcounts(x)[variable.genes,]})
mnn.out <- do.call(fastMNN, c(original, list(k=20, d=50, auto.merge=TRUE)))

# set column names
colnames(reducedDim(mnn.out)) = paste0("MNN_", 1:ncol(reducedDim(mnn.out)))

# add MNN correction results to Seurat object
panc.immune[["mnn"]] <- CreateDimReducObject(embeddings=reducedDim(mnn.out)[rownames(panc.immune@meta.data),], key="MNN_", assay=DefaultAssay(panc.immune))

# Run non-linear dimensional reduction (UMAP)
panc.immune %<>% RunUMAP(dims=1:50, reduction="mnn", n.components=3, seed.use=42, n.neighbors=30, n.epochs=200)

# Cluster the cells
panc.immune %<>% FindNeighbors(reduction="mnn", dims=1:50)
panc.immune %<>% FindClusters(resolution=seq(0.05,2,by=0.05), verbose=T)

# set cell identity
panc.immune %<>% SetIdent(value="RNA_snn_res.0.8")

# merge clusters
# C0 + C1 + C4                ==>  C0 (Macrophage)
# C2                          ==>  C1 (Macrophage eating islet cells)
# C3                          ==>  C2 (MAST)
# C5                          ==>  C3 (T-cells)
# C6                          ==>  C4 (B-cells)
merge.clust <- c(0,0,1,2,0,3,4)
names(merge.clust) <- 0:6

final.clust <- as.vector(merge.clust[as.vector(Idents(panc.immune))])
names(final.clust) <- names(Idents(panc.immune))
final.clust <- factor(final.clust, levels=0:4)

# add final clusters to meta data
panc.immune[["final.clust"]] <- final.clust

# set final clusters
Idents(panc.immune) <- final.clust

# save Seurat object
saveRDS(panc.immune, file.path(immune.infodir, 'panc.immune.rds'))

# -------- merge immune cell sub-clusters back to the main Seurat object -------- #

# rename immune cell sub-clusters, starting from 8
immune.clusters.renamed <- FetchData(panc.immune, vars='ident') %>% rownames_to_column('cellID') %>% 
  mutate(ident=as.character(as.numeric(as.vector(ident)) + 8))

# merge sub-clusters into the main Seurat object; rename endothelial cluster (8) to 7
combined.clusters <- FetchData(panc, vars='ident') %>% rownames_to_column('cellID') %>% 
  mutate(ident=as.vector(ident)) %>% dplyr::filter(! cellID %in% immune.clusters.renamed$cellID) %>% 
  mutate(ident=ifelse(ident == '8', '7', ident)) %>% 
  dplyr::bind_rows(immune.clusters.renamed) %>% column_to_rownames('cellID') %>% 
  dplyr::rename('combined.clust'='ident')

combined.clusters$combined.clust <- factor(combined.clusters$combined.clust, levels=0:12)

panc <- AddMetaData(panc, metadata=combined.clusters)

Idents(panc) <- 'combined.clust'

# C0 - beta cells
# C1: alpha cells
# C2: acinar cells
# C3: ductal cells
# C4: mesenchymal cells
# C5: delta cells
# C6: PP cells
# C7: endothelial cells
# C8: macrophages
# C9: DC cells
# C10: immune progenitor cells
# C11: T cells
# C12: B cells

saveRDS(panc, file.path(infodir,'panc.combined_labels.rds'))
