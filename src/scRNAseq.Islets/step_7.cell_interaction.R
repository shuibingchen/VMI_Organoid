# step_7.cell_interaction.R
# --- process scRNAseq data from islets samples ---
# step 7: conduct cell-cell interaction analysis
# Author: Tuo Zhang
# Date: 8/1/2024
# 

library(Seurat)
library(CellChat)
library(tidyverse)
library(patchwork)
library(RColorBrewer)
library(ComplexHeatmap)

options(stringsAsFactors = FALSE)

# folders
workdir <- "."
figdir <- file.path(workdir, "figure")
infodir <- file.path(workdir, "info")

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

# set the ligand-receptor interaction database
CellChatDB <- CellChatDB.human

# use all CellChatDB for cell-cell communication analysis
CellChatDB.use <- CellChatDB # simply use the default CellChatDB

# Load Seurat object (all cells, combined immune sub-cluster labels)
panc <- readRDS(file.path(infodir,'panc.combined_labels.rds'))

# prepare a metadata sheet
rename.clusters <- c('0'='Beta', '1'='Alpha', '2'='Acinar', '3'='Ductal', '4'='Mesenchymal', 
                     '5'='Delta', '6'='PP', '7'='Endo', '8'='Macrophage', '9'='DC',
                     '10'='Immune_progenitor', '11'='T', '12'='B')

meta <- FetchData(panc, vars=c('ident')) %>% mutate(ident=as.vector(ident)) %>% 
  mutate(ident=rename.clusters[ident]) %>% 
  dplyr::rename('labels'='ident')

# define cell groups
cov2.mock.samples <- c('control-9-mock','control-10-mock','control-11-mock')
cov2.cov2.samples <- c('control-9-SARS-CoV-2','control-10-SARS-CoV-2','control-11-SARS-CoV-2')
cvb4.mock.samples <- c('control-10-mock','control-12-mock')
cvb4.cvb4.samples <- c('control-10-CVB4','control-12-CVB4-1','control-12-CVB4-2')

# calculate number of cells per cluster (cell type) per group
stats.cov2.mock <- FetchData(panc, vars=c('Name','ident')) %>% dplyr::filter(Name %in% cov2.mock.samples) %>% 
  group_by(ident) %>% summarise_at('Name', list(nCells=length)) %>% mutate(group='MOCK_COV2', condition='MOCK') %>% 
  mutate(ident=as.numeric(as.character(ident))) %>% 
  tidyr::complete(ident=0:12, fill=list(nCells=0, group='MOCK_COV2', condition='MOCK'))

stats.cov2.cov2 <- FetchData(panc, vars=c('Name','ident')) %>% dplyr::filter(Name %in% cov2.cov2.samples) %>% 
  group_by(ident) %>% summarise_at('Name', list(nCells=length)) %>% mutate(group='MOCK_COV2', condition='COV2') %>% 
  mutate(ident=as.numeric(as.character(ident))) %>% 
  tidyr::complete(ident=0:12, fill=list(nCells=0, group='MOCK_COV2', condition='COV2'))

stats.cvb4.mock <- FetchData(panc, vars=c('Name','ident')) %>% dplyr::filter(Name %in% cvb4.mock.samples) %>% 
  group_by(ident) %>% summarise_at('Name', list(nCells=length)) %>% mutate(group='MOCK_CVB4', condition='MOCK') %>% 
  mutate(ident=as.numeric(as.character(ident))) %>% 
  tidyr::complete(ident=0:12, fill=list(nCells=0, group='MOCK_CVB4', condition='MOCK'))

stats.cvb4.cvb4 <- FetchData(panc, vars=c('Name','ident')) %>% dplyr::filter(Name %in% cvb4.cvb4.samples) %>% 
  group_by(ident) %>% summarise_at('Name', list(nCells=length)) %>% mutate(group='MOCK_CVB4', condition='CVB4') %>% 
  mutate(ident=as.numeric(as.character(ident))) %>% 
  tidyr::complete(ident=0:12, fill=list(nCells=0, group='MOCK_CVB4', condition='CVB4'))

stats <- rbind(stats.cov2.mock, stats.cov2.cov2, stats.cvb4.mock, stats.cvb4.cvb4)

# skip clusters with fewer than 50 cells
clusters.to.skip <- as.character(stats %>% dplyr::filter(nCells <= 50) %>% pull(ident) %>% unique())
clusters.use <- setdiff(as.character(unique(Idents(panc))), clusters.to.skip)

# Run CellChat
cellchat.cov2.mock <- my.preprocessing.cellchat(seurat.obj=panc, sample.use=cov2.mock.samples, clusters.use=clusters.use, 
                                                metadata=meta, group.by='labels', CellChatDB.use=CellChatDB.use, 
                                                population.size=TRUE, infodir=infodir, 
                                                suffix='cov2.mock', sources.use = c(3,7), targets.use = c(3,7), workers=4)

cellchat.cov2.cov2 <- my.preprocessing.cellchat(seurat.obj=panc, sample.use=cov2.cov2.samples, clusters.use=clusters.use, 
                                                metadata=meta, group.by='labels', CellChatDB.use=CellChatDB.use, 
                                                population.size=TRUE, infodir=infodir, 
                                                suffix='cov2.cov2', sources.use = c(3,7), targets.use = c(3,7), workers=4)

cellchat.cvb4.mock <- my.preprocessing.cellchat(seurat.obj=panc, sample.use=cvb4.mock.samples, clusters.use=clusters.use, 
                                                metadata=meta, group.by='labels', CellChatDB.use=CellChatDB.use, 
                                                population.size=TRUE, infodir=infodir, 
                                                suffix='cvb4.mock', sources.use = c(3,7), targets.use = c(3,7), workers=4)

cellchat.cvb4.cvb4 <- my.preprocessing.cellchat(seurat.obj=panc, sample.use=cvb4.cvb4.samples, clusters.use=clusters.use, 
                                                metadata=meta, group.by='labels', CellChatDB.use=CellChatDB.use, 
                                                population.size=TRUE, infodir=infodir, 
                                                suffix='cvb4.cvb4', sources.use = c(3,7), targets.use = c(3,7), workers=4)

# merge MOCK and COV2 cellchat objects
object.list.cov2 <- list(MOCK = cellchat.cov2.mock, COV2 = cellchat.cov2.cov2)
cellchat.cov2.merged <- mergeCellChat(object.list.cov2, add.names = names(object.list.cov2))

# merge MOCK and CVB4 cellchat objects
object.list.cvb4 <- list(MOCK = cellchat.cvb4.mock, CVB4 = cellchat.cvb4.cvb4)
cellchat.cvb4.merged <- mergeCellChat(object.list.cvb4, add.names = names(object.list.cvb4))

# save merged CellChat objects
saveRDS(cellchat.cov2.merged, file.path(infodir, "cellchat.cov2.merged.rds"))
saveRDS(cellchat.cvb4.merged, file.path(infodir, "cellchat.cvb4.merged.rds"))
