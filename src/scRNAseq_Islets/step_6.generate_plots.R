# step_6.generate_plots.R
# --- process scRNAseq data from islets samples ---
# step 6: generate plots in the manuscript
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
solo.infodir <- file.path(infodir, 'solo')

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

# load Seurat object (all cells)
panc <- readRDS(file.path(infodir, "panc.rds"))

# set color
my.cluster.color <- c(brewer.pal(12,'Paired')[1:9])
names(my.cluster.color) <- 0:8

# set x/y-axis boundaries
x.upper <- 11
x.lower <- -7
y.upper <- 6
y.lower <- -14

# Figure S2A: UMAP of human islets exposed to mock, SARS-CoV-2 or CVB4 viruses
g <- myDimPlot(tobj=panc, treduct="umap", tcate="ident", tsuffix="Cluster", tcolor=my.cluster.color, tlabel=FALSE, tsplit=FALSE, 
               txlim=c(x.lower,x.upper), tylim=c(y.lower,y.upper), tptsize=1.5, talpha=0.7, tltsize=20, tatlsize=22)
ggsave(file.path(figdir, "Fig.S2A.png"), plot=g, width=10, height=8, dpi=300)

# Figure S2B: Violin plot of cell markers of each cell population
known.markers <- c('INS','GCG','PRSS2','KRT19','COL1A1','SST','PPY','IFI30','PECAM1')
for (tgene in known.markers){
  print(tgene)
  plot <- MyExpViolin(tobj=panc, tgene=tgene, tgroup_by='ident', tgroup_order=0:8, tcolor_by='ident', tcolor_order=0:8, 
  	tcolor=my.cluster.color, tcells=NULL, tassay="RNA", tncol=1)
  ggsave(file.path(figdir, paste("Fig.S2B", tgene, "png", sep='.')), plot=plot, height=4, width=6, dpi=300)
}

# Figure S2C: UMAP showed the expression of CVB4 polyprotein in human islets exposed to mock or CVB4 virus
cvb4.genes <- grep(cvb4.pattern, rownames(panc), value=T)

# subset CVB4 mock cells
cvb4.mock.cells <- FetchData(panc, vars=c('Name')) %>% rownames_to_column('cellID') %>% 
  dplyr::filter(Name %in% c("control-10-mock","control-12-mock")) %>% pull('cellID')
panc.cvb4.mock <- subset(panc, cells=cvb4.mock.cells)

g <- FeaturePlot(panc.cvb4.mock, cvb4.genes)
ggsave(file.path(figdir, "Fig.S2C-1.png"), plot=g, width=8, height=6.5, dpi=300)

# subset CVB4 infected cells
cvb4.infected.cells <- FetchData(panc, vars=c('Name')) %>% rownames_to_column('cellID') %>% 
  dplyr::filter(Name %in% c('control-10-CVB4',"control-12-CVB4-1","control-12-CVB4-2")) %>% pull('cellID')
panc.cvb4.infect <- subset(panc, cells=cvb4.infected.cells)

g <- FeaturePlot(panc.cvb4.infect, cvb4.genes)
ggsave(file.path(figdir, "Fig.S2C-2.png"), plot=g, width=8, height=6.5, dpi=300)

# Figure S2D: Jitter plot showed the expression of CVB4 polyprotein in human islets exposed to mock or CVB4 virus
tdata <- FetchData(panc.cvb4.infect, vars=c('ident',cvb4.genes))
g <- ggplot(tdata, aes(x=ident, y=polyprotein))
g <- g + geom_jitter(aes(color=ident), position=position_jitter(0.2), size=1, alpha=0.8)
g <- g + scale_color_manual(values=my.cluster.color)
g <- g + theme_bw()
g <- g + theme(axis.title=element_text(size=15), axis.text=element_text(size=13))
g <- g + ylab("Expression Level") + xlab("Identity")
ggsave(file.path(figdir, "Fig.S2D.png"), plot=g, width=7, height=4, dpi=300)

# load Seurat object (immune cells)
panc.immune <- readRDS(file.path(infodir, "panc.immune.rds"))

# Figure 2A: UMAP of immune cell populations in human islets exposed to mock, SARS-CoV-2 or CVB4

# set x/y-axis boundaries
x.upper <- 3.5
x.lower <- -12
y.upper <- 6
y.lower <- -3

# set color
my.immune.cluster.color <- c(brewer.pal(12,'Paired')[c(2,4,12,10,6)])
names(my.immune.cluster.color) <- 0:4

g <- myDimPlot(tobj=panc.immune, treduct="umap", tcate="ident", tsuffix="Cluster", tcolor=my.immune.cluster.color, tlabel=TRUE, tsplit=FALSE, 
               txlim=c(x.lower,x.upper), tylim=c(y.lower,y.upper), tptsize=1.5, talpha=0.7, tltsize=20, tatlsize=22)
ggsave(file.path(figdir, "Fig.2A.png"), plot=g, width=10, height=8, dpi=300)

# Figure 2B: UMAP and violin plots of immune cell markers

# set color for MyFeaturePlot
myExpLowColor <- '#d9d9d9'
myExpHighColor <- '#b30000'

# markers
immune.markers <- c('CD14','CD3E','CD19','IL3RA','KIT')

for (tgene in immune.markers){
  plot <- MyFeaturePlot(tobj=panc.immune, tgenes=c(tgene), tcells=NULL, tassay="RNA", treduction.name="umap", tlowcolor=myExpLowColor, thighcolor=myExpHighColor,
                        txlim=c(x.lower,x.upper), tylim=c(y.lower,y.upper), tncol=1, tlegend=NULL)
  ggsave(file.path(figdir, "Fig.2B-1.png"), plot=plot, width=10, height=8, dpi=300)
  plot <- MyExpViolin(tobj=panc.immune, tgene=tgene, tgroup_by='ident', tgroup_order=0:4, tcolor_by='ident', tcolor_order=0:4, 
                      tcolor=my.immune.cluster.color, tcells=NULL, tassay="RNA", tncol=1)
  ggsave(file.path(figdir, "Fig.2B-2.png"), plot=plot, height=4, width=6, dpi=300)
}
