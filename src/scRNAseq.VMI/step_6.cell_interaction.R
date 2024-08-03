# step_6.cell_interaction.R
# --- process scRNAseq data from VMI samples ---
# step 6: Conduct cell-cell interactions analysis between macrophage and beta cells
# Author: Tuo Zhang
# Date: 8/1/2024
# 

library(Seurat)
library(CellChat)
library(tidyverse)
library(patchwork)
library(magrittr)
library(future)
library(RColorBrewer)
library(ComplexHeatmap)

# folders
workdir <- "."
sourcedir <- file.path(workdir, "data")
figdir <- file.path(workdir, "figure")
infodir <- file.path(workdir, "info")

# random seed
rseed <- 98
set.seed(rseed)

options(stringsAsFactors = FALSE)

# set parallelization in Seurat
plan("multicore", workers=6)
options(future.globals.maxSize=10*1024^3)

# set the ligand-receptor interaction database
CellChatDB <- CellChatDB.human

# use all CellChatDB for cell-cell communication analysis
CellChatDB.use <- CellChatDB # simply use the default CellChatDB

# load RNA Seurat object
panc <- readRDS(file.path(rna.infodir, 'panc.rds'))

# Subset seurat object for beta and macrophage cells only
panc.bm <- subset(panc, idents=c(0,8))

# Rename identity classes
panc.bm <- RenameIdents(panc.bm, '0'='0', '1'='8')

# prepare a metadata sheet
rename.clusters <- c('0'='Beta', '1'='Macrophage')

meta <- FetchData(panc.bm, vars=c('ident')) %>% mutate(ident=as.vector(ident)) %>% 
  mutate(ident=rename.clusters[ident]) %>% 
  dplyr::rename('labels'='ident')

# preprocess BUC
cellchat.BUC <- my.preprocessing.cellchat(seurat.obj=panc.bm, sample.use=c('BUC'), 
                                          clusters.use=c(0,1), metadata=meta, group.by='labels',
                                          CellChatDB.use=CellChatDB.use, 
                                          population.size=TRUE, infodir=infodir.1, 
                                          suffix='bm.BUC', 
                                          sources.use = c(1,2), targets.use = c(1,2), workers=4)

# preprocess BM0
cellchat.BM0 <- my.preprocessing.cellchat(seurat.obj=panc.bm, sample.use=c('BM0'), 
                                          clusters.use=c(0,1), metadata=meta, group.by='labels',
                                          CellChatDB.use=CellChatDB.use, 
                                          population.size=TRUE, infodir=infodir.1, 
                                          suffix='bm.BM0', 
                                          sources.use = c(1,2), targets.use = c(1,2), workers=4)

# preproces BM1
cellchat.BM1 <- my.preprocessing.cellchat(seurat.obj=panc.bm, sample.use=c('BM1'), 
                                          clusters.use=c(0,1), metadata=meta, group.by='labels',
                                          CellChatDB.use=CellChatDB.use, 
                                          population.size=TRUE, infodir=infodir.1, 
                                          suffix='bm.BM1', 
                                          sources.use = c(1,2), targets.use = c(1,2), workers=4)

# merge cellchat objects
object.list <- list(BUC = cellchat.BUC, BM0 = cellchat.BM0, BM1 = cellchat.BM1)
cellchat.merged <- mergeCellChat(object.list, add.names = names(object.list))

# save object
saveRDS(cellchat.merged, file.path(infodir, 'cellchat.merged.rds'))
