# step_6.generate_plots.R
# --- process scRNAseq data from islets samples ---
# step 6: Generate plots in the manuscript
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

# Figure 4B: Dot plot displaying cell markers of each cluster using scRNA-seq dataset
cell.type.markers <- rev(c('INS','SST','NEUROG3','GCG','KRT19','PRSS1','PECAM1','CD14'))

plot <- DotPlot(panc, features=cell.type.markers) + coord_flip()
ggsave(file.path(figdir, 'Fig.4B.png'), width=7, height=5, plot=plot, dpi=300)

# Figure 4E: Volcano plot of DE genes in β cell cluster of VMI organoids at day 7 after reaggregation versus separately cultured cells as analyzed by scRNA-seq
plot <- myVolcanoPlot(tdefile=file.path(infodir, 'DE.C0-BM0.vs.C0-BUC.mincells_10.wilcox.min_pct_0.1.logfc_0.1.txt'), 
	tx='p_val', ty='avg_logFC', tcutFC=0.25, tcutP=20, tlabel.genes.up=c('RPL13A','SMIM32'), tlabel.genes.down=c('AFP','GCG','ACTB','SST','PRSS2'), 
	txlabel='LogFoldChange', tylabel='-LogP-value', tupper=300, talpha=0.8, tcolor.up='firebrick3', tcolor.down='steelblue3', tcolor.other='gray60')
ggsave(file.path(figdir, 'Fig.4F.png'), width=7, height=5, plot=plot, dpi=300)

