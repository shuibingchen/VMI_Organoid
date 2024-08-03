# step_7.generate_plots.R
# --- process scRNAseq data from VMI samples ---
# step 7: Generate plots in the manuscript
# Author: Tuo Zhang
# Date: 8/1/2024
# 

library(Seurat)
library(CellChat)
library(tidyverse)
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

# sample info
sample.info <- data.frame(SeqName=c('non-co-RNA','VMI-UN-RNA','VMI-PRO-RNA'),
                          Group=c('non-coculture','VMI_unstimulated','VMI_pro-inflammatory'),
                          Name=c('BUC','BM0','BM1'))
rownames(sample.info) <- sample.info$Name

# Load Seurat object
panc <- readRDS(file.path(infodir, 'pand.rds'))

# Load cellchat object
cellchat.merged <- readRDS(file.path(infodir, 'cellchat.merged.rds'))

# Figure 4B: Dot plot displaying cell markers of each cluster using scRNA-seq dataset
cell.type.markers <- rev(c('INS','SST','NEUROG3','GCG','KRT19','PRSS1','PECAM1','CD14'))

plot <- DotPlot(panc, features=cell.type.markers) + coord_flip()
ggsave(file.path(figdir, 'Fig.4B.png'), width=7, height=5, plot=plot, dpi=300)

# Figure 4E: Volcano plot of DE genes in β cell cluster of VMI organoids at day 7 after reaggregation versus separately cultured cells as analyzed by scRNA-seq
plot <- myVolcanoPlot(tdefile=file.path(infodir, 'DE.C0-BM0.vs.C0-BUC.mincells_10.wilcox.min_pct_0.1.logfc_0.1.txt'), 
	tx='p_val', ty='avg_logFC', tcutFC=0.25, tcutP=20, tlabel.genes.up=c('RPL13A','SMIM32'), tlabel.genes.down=c('AFP','GCG','ACTB','SST','PRSS2'), 
	txlabel='LogFoldChange', tylabel='-LogP-value', tupper=300, talpha=0.8, tcolor.up='firebrick3', tcolor.down='steelblue3', tcolor.other='gray60')
ggsave(file.path(figdir, 'Fig.4F.png'), width=7, height=5, plot=plot, dpi=300)

# Figure 4F: Dot plot analysis of β cell associated genes in β cell cluster of VMI organoids at day 7 after reaggregation and separately cultured cells as analyzed by scRNA-seq.
gene.set.1 <- c("INS","PDX1","PAX6","HNF1B","PIK3CB","SLC2A1")

# extract cell IDs for beta cells from BM0 and BUC samples
beta.BM0.BUC.cells <- rownames(FetchData(panc, vars=c('merged.cluster','Name')) %>% 
                               dplyr::filter(merged.cluster == 0 & Name %in% c('BUC','BM0')))

panc.beta.BM0.BUC <- subset(panc, cells=beta.BM0.BUC.cells)
panc.beta.BM0.BUC$Name <- factor(panc.beta.BM0.BUC$Name, levels=c('BUC','BM0'))

g <- DotPlot(panc.beta.BM0.BUC, features=gene.set.1, group.by='Name') + coord_flip()
ggsave(file.path(figdir, 'Fig.4F.png'), width=5, height=5, dpi=300)

# Figure 4H: Dot plot analysis of endothelial cell associated genes in endothelial cell cluster of VMI organoids at day 7 after reaggregation and separately cultured cells as analyzed by scRNA-seq.
gene.set.2 <- c('INSR','VWF','PDGFB','EDN1','S1PR1','RSPO3')

# extract cell IDs for endothelial cells from BM0 and BUC samples
endo.BM0.BUC.cells <- rownames(FetchData(panc, vars=c('merged.cluster','Name')) %>% 
                               dplyr::filter(merge.clust.1 == 7 & Name %in% c('BM0','BUC')))

g <- DotPlot(subset(panc, cells=endo.BM0.BUC.cells), features=gene.set.2, group.by='Name') + coord_flip()
ggsave(file.path(figdir, 'Fig.4H.png'), width=5, height=4, dpi=300)

# Figure 6A: Dot plot showed the differential signaling from macrophages to β cells in VMI organoids containing unstimulated or proinflammatory macrophages at day 7 after reaggregation.
png(file.path(figdir, 'Fig.6A.png'), width=4, height=4.5, units='in', res=300)
netVisual_bubble(cellchat.merged, sources.use=c(2), targets.use = c(1),  
                 comparison = c(2, 3), angle.x = 45, line.on=F)
dev.off()
