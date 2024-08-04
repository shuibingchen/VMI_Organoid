# step_5.generate_plots.R
# --- process scRNAseq data from VMI samples ---
# step 5: Generate plots in the manuscript
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
panc.atac <- readRDS(file.path(infodir, 'panc.integrated.rds'))

# load coembeded object
coembed.v4 <- readRDS(file.path(infodir, 'coembed.v4.rds'))

# prepare RNA cell labels
tlabel.rna <- FetchData(panc.rna, vars=c('merged.cluster','Name')) %>% rownames_to_column('cellID') %>% 
    mutate(cellID=paste0('RNA_', cellID)) %>% column_to_rownames('cellID') %>% 
    mutate(transfer.cluster=paste0('C',merged.cluster), 
           Name=paste0('RNA_',Name), dataset='RNA')

# prepare ATAC cell labels
tlabel.atac <- FetchData(panc.atac, vars=c('dataset','predicted.id')) %>% 
    rownames_to_column('cellID') %>% mutate(cellID=paste0('ATAC_', cellID)) %>% 
    column_to_rownames('cellID') %>% 
    rename(c('transfer.cluster'='predicted.id','Name'='dataset')) %>% 
    mutate(transfer.cluster=paste0('C',transfer.cluster), 
           Name=paste0('ATAC_',Name), dataset='ATAC')

# combine cell labels and add to the coembeded object
coembed.v4 %<>% AddMetaData(metadata=rbind(tlabel.rna, tlabel.atac))

# set color
my.cluster.color <- brewer.pal(12,'Paired')[1:9]
names(my.cluster.color) <- paste0('C',0:8)

# Figure 4A: Integrative UMAP of scRNA-seq and snATAC-seq analysis of VMI organoids 
plot <- myDimPlot3(tobj=coembed.v4, treduct="umap", 
                   tgroup_by="transfer.cluster", tgroup_order=paste0('C',0:8), tsuffix="Cluster", 
                   tcolor=my.cluster.color, tsplit_by='dataset', tsplit_order=c('ATAC','RNA'), 
                   tlabel=T, tncol=2, tptsize=0.4, tlbsize=3.5) + theme(legend.position="none")
ggsave(file.path(figdir, "Fig.4A.png"), 
       plot=plot, width=8.5, height=4, dpi=300)

# Figure 4C: Individual UMAP of scRNA-seq and snATAC-seq analysis of VMI organoids
cells.BUC.BM0 <- FetchData(coembed.v4, vars=c('Name')) %>% rownames_to_column('cellID') %>%
    dplyr::filter(Name %in% c('RNA_BUC','RNA_BM0','ATAC_BUC','ATAC_BM0')) %>% 
    pull('cellID')

plot <- myDimPlot3(tobj=subset(coembed.v4, cells=cells.BUC.BMO), treduct="umap", 
                   tgroup_by="transfer.cluster", tgroup_order=paste0('C',0:8), tsuffix="Cluster", 
                   tcolor=my.cluster.color, tsplit_by='Name', 
                   tsplit_order=c(paste('ATAC', c('BUC','BM0'), sep='_'), 
                                  paste('RNA', c('BUC','BM0'), sep='_')), 
                   tlabel=T, tncol=3, tptsize=0.4, tlbsize=3.5) + theme(legend.position="none")
ggsave(file.path(figdir, "Fig.4C.png"), plot=plot, width=8.5, height=8, dpi=300)

# Figure 4D: Pie chart showed the relative percentages of each cell types in VMI organoids
stat.RNA_BM0 <- FetchData(coembed.v4, vars=c('Name','transfer.cluster')) %>% 
    dplyr::filter(Name == 'RNA_BM0') %>% group_by(transfer.cluster) %>% count() %>% as.data.frame() %>% 
    tidyr::complete(transfer.cluster=factor(paste0('C',0:8)), fill=list(n=0)) %>%
    mutate_at('n', list(percent=function(x) { round(x/sum(x)*100,2) }))

write.table(stat.RNA_BM0, file.path(infodir, "Fig.4D.RNA.txt"), quote=F, sep='\t')

stat.ATAC_BM0 <- FetchData(coembed.v4, vars=c('Name','transfer.cluster')) %>% 
    dplyr::filter(Name == 'ATAC_BM0') %>% group_by(transfer.cluster) %>% count() %>% as.data.frame() %>% 
    tidyr::complete(transfer.cluster=factor(paste0('C',0:8)), fill=list(n=0)) %>%
    mutate_at('n', list(percent=function(x) { round(x/sum(x)*100,2) }))

write.table(stat.ATAC_BM0, file.path(infodir, "Fig.4D.ATAC.txt"), quote=F, sep='\t')

blank_theme <- theme_minimal()+
  theme(
  axis.title.x = element_blank(),
  axis.title.y = element_blank(),
  panel.border = element_blank(),
  panel.grid=element_blank(),
  axis.ticks = element_blank(),
  plot.title=element_text(size=14, face="bold")
  )

p1 <- ggplot(stat.RNA_BM0, aes(x="", y=percent, fill=transfer.cluster))
p1 <- p1 + geom_bar(width = 1, stat="identity") + coord_polar("y")
p1 <- p1 + scale_fill_manual(values=my.cluster.color) +  blank_theme
p1 <- p1 + ggtitle('BM0 - RNA')

p2 <- ggplot(stat.ATAC_BM0, aes(x="", y=percent, fill=transfer.cluster))
p2 <- p2 + geom_bar(width = 1, stat="identity") + coord_polar("y")
p2 <- p2 + scale_fill_manual(values=my.cluster.color) +  blank_theme
p2 <- p2 + ggtitle('BM0 - ATAC')

g <- p1 | p2
ggsave(file.path(figdir, 'Fig.4D.png'), plot=g, width=9.5, height=4.5, dpi=300)

# Figure 4G: Chromatin accessibility signals of SLC2A1, INS, PDX1 in the β cell cluster of VMI organoids at day 7 after reaggregation and separately cultured cells as analyzed by snATAC-seq. 
# 1. SLC2A1
plot <- my.coverage.plot.beta.cmp(tgene='SLC2A1', cdt=c('BUC','BM0'), tsuffix='BM0_vs_BUC', 
    up=2000, down=2000, include.tile=T)
ggsave(file.path(figdir, 'Fig.4D.SLC2A1.png'), plot=g, width=5, height=5, dpi=300)
# 2. INS
plot <- my.coverage.plot.beta.cmp(tgene='INS', cdt=c('BUC','BM0'), tsuffix='BM0_vs_BUC', 
    up=25000, down=25000, include.tile=T)
ggsave(file.path(figdir, 'Fig.4D.INS.png'), plot=g, width=5, height=5, dpi=300)
# 3. PDX1
plot <- my.coverage.plot.beta.cmp(tgene='PDX1', cdt=c('BUC','BM0'), tsuffix='BM0_vs_BUC', 
    up=15000, down=5000, include.tile=T)
ggsave(file.path(figdir, 'Fig.4D.PDX1.png'), plot=g, width=5, height=5, dpi=300)

# Figure 5D: Integrative UMAP of VMI organoids at day 7 after reaggregation containing unstimulated or proinflammatory macrophages.
cells.BM0.BM1 <- FetchData(coembed.v4, vars=c('Name')) %>% rownames_to_column('cellID') %>%
    dplyr::filter(Name %in% c('RNA_BM0','RNA_BM1','ATAC_BM0','ATAC_BM1')) %>% 
    pull('cellID')

plot <- myDimPlot3(tobj=subset(coembed.v4, cells=cells.BM0.BM1), treduct="umap", 
                   tgroup_by="transfer.cluster", tgroup_order=paste0('C',0:8), tsuffix="Cluster", 
                   tcolor=my.cluster.color, tsplit_by='Name', 
                   tsplit_order=c(paste('ATAC', c('BM0','BM1'), sep='_'), 
                                  paste('RNA', c('BM0','BM1'), sep='_')), 
                   tlabel=T, tncol=3, tptsize=0.4, tlbsize=3.5) + theme(legend.position="none")
ggsave(file.path(figdir, "Fig.5D.png"), plot=plot, width=8.5, height=8, dpi=300)

# Figure 5I: Chromatin accessibility signals of CASP1, CASP9, IL1B and NLRP3 in the β cell cluster of VMI organoids at day 7 after reaggregation containing unstimulated or proinflammatory macrophages.
# 1. CASP1
plot <- my.coverage.plot.beta.cmp(tgene='CASP1', cdt=c('BM0','BM1'), tsuffix='BM1_vs_BM0', 
    up=2000, down=2000, width=5, height=2.5, include.tile=T)
ggsave(file.path(figdir, 'Fig.5I.CASP1.png'), plot=g, width=5, height=2.5, dpi=300)
# CASP9
plot <- my.coverage.plot.beta.cmp(tgene='CASP9', cdt=c('BM0','BM1'), tsuffix='BM1_vs_BM0', 
    up=-10000, down=-15000, width=5, height=2.5, include.tile=T)
ggsave(file.path(figdir, 'Fig.5I.CASP9.png'), plot=g, width=5, height=2.5, dpi=300)
# IL1B
plot <- my.coverage.plot.beta.cmp(tgene='IL1B', cdt=c('BM0','BM1'), tsuffix='BM1_vs_BM0', 
    up=2000, down=5000, width=5, height=2.5, include.tile=T)
ggsave(file.path(figdir, 'Fig.5I.IL1B.png'), plot=g, width=5, height=2.5, dpi=300)
# NLRP3
plot <- my.coverage.plot.beta.cmp(tgene='NLRP3', cdt=c('BM0','BM1'), tsuffix='BM1_vs_BM0', 
    up=10000, down=2000, width=5, height=2.5, include.tile=T)
ggsave(file.path(figdir, 'Fig.5I.NLRP3.png'), plot=g, width=5, height=2.5, dpi=300)

