# step_2.doublets_detection.R
# --- process scRNAseq data from VMI samples ---
# step 2: Identify doublets
# Author: Tuo Zhang
# Date: 8/1/2024
# 

library(DoubletFinder)
library(Seurat)
library(tidyverse)
library(Matrix)
library(magrittr)

# folders
workdir <- "."
sourcedir <- file.path(workdir, "data")
figdir <- file.path(workdir, "figure")
infodir <- file.path(workdir, "info")

# Set a random seed
set.seed(98)

# Load Seurat object
panc <- readRDS(file.path(infodir, 'panc.initial.rds'))

# doublet rate
# https://assets.ctfassets.net/an68im79xiti/1eX2FPdpeCgnCJtw4fj9Hx/7cb84edaa9eca04b607f9193162994de/CG000204_ChromiumNextGEMSingleCell3_v3.1_Rev_D.pdf
# according to 10X Genomics document, the doublet rate is correlated with the number of targeted/capture cells
# 1000 recovered cells ~ 0.8%
doublet.rate <- 0.008

# number of PCs
pcs <- 15

# resolution
res <- 0.8

# pN: number of generated artificial doublets, as a proportion of the merged real-artificial data
pN <- 0.25

# doublet results
df.results <- data.frame(cellID=character(), orig.ident=character(), df.isdoublet=character(), df.score=numeric())

# detect doublet per sample
for (sid in unique(panc$orig.ident)){
    # subset Seurat object
    seu <- subset(panc, subset=orig.ident %in% c(sid))
    # pre-process Seurat object
    seu %<>% NormalizeData()
    seu %<>% FindVariableFeatures(selection.method="vst", nfeatures=2000)
    seu %<>% ScaleData()
    seu %<>% RunPCA()
    seu %<>% RunUMAP(dims=1:pcs)
    seu %<>% FindNeighbors(reduction="pca", dims=1:pcs)
    seu %<>% FindClusters(resolution=res)
    # pK Identification (no ground-truth)
    sweep.res.list <- paramSweep_v3(seu, PCs=1:pcs, sct=FALSE)
    sweep.stats <- summarizeSweep(sweep.res.list, GT=FALSE)
    bcmvn <- find.pK(sweep.stats)
    # factor --> numeric
    bcmvn$pK <- as.numeric(as.character(bcmvn$pK))
    # select optimal pK
    pK <- bcmvn$pK[which(bcmvn$BCmetric %in% max(bcmvn$BCmetric))]
    # Homotypic Doublet Proportion Estimate
    homotypic.prop <- modelHomotypic(seu@meta.data$seurat_clusters)
    nExp_poi <- round(doublet.rate*ncol(seu)/1000*ncol(seu))
    nExp_poi.adj <- round(nExp_poi*(1-homotypic.prop))
    # run doublet detection
    seu %<>% doubletFinder_v3(PCs=1:pcs, pN=pN, pK=pK, nExp=nExp_poi, reuse.pANN=FALSE, sct=FALSE)
    seu %<>% doubletFinder_v3(PCs=1:pcs, pN=pN, pK=pK, nExp=nExp_poi.adj, reuse.pANN=paste("pANN",pN,pK,nExp_poi,sep="_"), sct=FALSE)
    # save pNN, singlet/doublet classification (original+adjusted), seurat_clusters (RNA_snn_res.0.8), UMAP coordinates(new)
    tdata <- FetchData(seu, vars=c(paste("pANN",pN,pK,nExp_poi,sep="_"),
                                   paste("DF.classifications",pN,pK,nExp_poi,sep="_"),
                                   paste("DF.classifications",pN,pK,nExp_poi.adj,sep="_")))
    colnames(tdata) <- c('df.score','DF.withHomo','df.isdoublet')
    tdata <- tdata %>% rownames_to_column('cellID') %>% mutate(orig.ident=sid) %>% 
    	dplyr::select('cellID', 'orig.ident', 'df.isdoublet', 'df.score')
    # add to results table
    df.results <- rbind(df.results, tdata)
}

write.table(df.results, file=file.path(infodir, "doublet.results.txt"), quote=F, sep='\t', row.names=F)
