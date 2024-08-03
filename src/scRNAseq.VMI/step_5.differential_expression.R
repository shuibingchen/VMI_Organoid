# step_5.differential_expression.R
# --- process scRNAseq data from islets samples ---
# step 5: perform differential expression analysis
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

# Load Seurat object
panc <- readRDS(file.path(infodir, 'pand.rds'))

# define comparison groups
tinfo <- FetchData(panc, vars=c('ident','Name')) %>% 
    mutate(compare=paste0('C', ident, '-', Name))

# add to the seurat object
panc %<>% AddMetaData(metadata=tinfo[,'compare',drop=F])

# set compare clusters
Idents(panc) <- 'compare'

# all possible groups
all.groups <- c()
for (clust in 0:8){
  all.groups <- c(all.groups, paste0('C',clust,'-',c('BUC','BM0','BM1')))
}

# get number of cells per group
nCells <- FetchData(panc, vars=c('ident')) %>% dplyr::count(ident) %>% tidyr::complete(ident=all.groups, fill=list(n=0))

# for each cluster, set up DE analysis
for (clust in 0:8){
    # BM0 v.s BUC
    my.DE.pair(panc, paste0('C',clust,'-BM0'), paste0('C',clust,'-BUC'), ntop=20, nCells=nCells,
               ttitle=paste0('C',clust,' co-culture M0 vs. un-coculture'), 
               no.mito=FALSE, mincells=10, min.pct=0.1, logfc.threshold=0.1, test.use='wilcox',
               cut.padj=0.1, cut.avglogfc=0, cutFC=0.6, cutP=20, 
               tfigdir=figdir, tinfodir=infodir)
    # BM1 v.s BM0
    my.DE.pair(panc, paste0('C',clust,'-BM1'), paste0('C',clust,'-BM0'), ntop=20, nCells=nCells,
               ttitle=paste0('C',clust,' co-culture M1 vs. co-culture M0'), 
               no.mito=FALSE, mincells=10, min.pct=0.1, logfc.threshold=0.1, test.use='wilcox',
               cut.padj=0.1, cut.avglogfc=0, cutFC=0.6, cutP=20, 
               tfigdir=figdir, tinfodir=infodir)
}

# set original clusters back
Idents(panc) <- 'merged.cluster'
