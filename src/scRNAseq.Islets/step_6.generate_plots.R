# step_8.generate_plots.R
# --- process scRNAseq data from islets samples ---
# step 8: generate plots in the manuscript
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
library(CellChat)

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

# Figure 2C: Dot plot analysis of proinflammatory macrophage-associated genes in macrophages of human islets exposed to mock or SARS-CoV-2 (MOI=1). 
m1.macrophage.markers <- c('IL1B','IL6','CXCL8','TNF','CCL2','CXCL10','IDO1','CD80')

macrophage.Mock.COV2.cells <- rownames(subset(FetchData(panc.immune, vars=c('ident','Name')), 
                                              ident %in% c(0) & Name %in% c('control-9-mock','control-9-SARS-CoV-2','control-10-mock','control-10-SARS-CoV-2',
                                                'control-11-mock','control-11-SARS-CoV-2')))
plot <- DotPlot.2(subset(panc.immune, cells=macrophage.Mock.COV2.cells), assay='RNA', features=m1.macrophage.markers, 
                  cols="RdBu", group.by='Condition', order.ids=c('Mock','COV2')) + 
  theme_bw() + coord_flip() + theme(axis.text.x=element_text(angle=60, hjust=1))
ggsave(file.path(figdir, "Fig.2C.png"), plot=plot, width=3.5, height=4.5, dpi=300)

# Figure 2G: Dot plot analysis of pyroptosis associated genes in the β cell cluster of human islets exposed to mock or SARS-CoV-2 (MOI=1).
gene.set.cov2 <- c("GSDME","GSDMD","CASP1","CASP8","NLRP3","IL18","CASP9","PYCARD","STS")

beta.cells <- rownames(subset(FetchData(panc, vars=c('ident')), ident %in% c(0)))

cov2.cells <- FetchData(panc, vars=c('Name')) %>% rownames_to_column('cellID') %>% 
  dplyr::filter(Name %in% c('control-9-mock','control-9-SARS-CoV-2','control-10-mock','control-10-SARS-CoV-2','control-11-mock','control-11-SARS-CoV-2')) %>% pull('cellID')

plot <- DotPlot.2(subset(panc, cells=intersect(beta.cells, cov2.cells)), assay='RNA', features=gene.set.cov2, cols="RdBu", group.by='Condition', 
                  order.ids=c('Mock','COV2')) + 
  theme_bw() + coord_flip() + theme(axis.text.x=element_text(angle=60, hjust=1))
ggsave(file.path(figdir, "Fig.2G.png"), plot=plot, width=3.5, height=4.5, dpi=300)

# Figure 2J: Dot plot analysis of proinflammatory macrophage-associated genes in the macrophage cluster of human islets exposed to mock or CVB4 
macrophage.Mock.CVB4.cells <- rownames(subset(FetchData(panc.immune, vars=c('ident','Name')), 
                                              ident %in% c(0) & Name %in% c('control-10-mock','control-10-CVB4','control-12-mock','control-12-CVB4-1','control-12-CVB4-2')))

plot <- DotPlot.2(subset(panc.immune, cells=macrophage.Mock.CVB4.cells), assay='RNA', features=m1.macrophage.markers, 
                  cols="RdBu", group.by='Condition', order.ids=c('Mock','CVB4')) + 
  theme_bw() + coord_flip() + theme(axis.text.x=element_text(angle=60, hjust=1))
ggsave(file.path(figdir, "Fig.2J.png",sep="."), plot=plot, width=3.5, height=4.5, dpi=300)

# Figure 2M: Dot plot analysis of pyroptosis pathway associated genes in the β cell cluster of human islets exposed to mock or CVB4 (2x106 PFU/ml).
gene.set.cvb4 <- c("GSDME","CASP1","NLRP3","IL1B","CASP9")

cvb4.cells <- FetchData(panc, vars=c('Name')) %>% rownames_to_column('cellID') %>% 
  dplyr::filter(Name %in% c('control-10-mock','control-10-CVB4','control-12-mock','control-12-CVB4-1','control-12-CVB4-2')) %>% pull('cellID')

plot <- DotPlot.2(subset(panc, cells=intersect(beta.cells, cvb4.cells)), assay='RNA', features=gene.set.cvb4, cols="RdBu", group.by='Condition', 
                  order.ids=c('Mock','CVB4')) + 
  theme_bw() + coord_flip() + theme(axis.text.x=element_text(angle=60, hjust=1))
ggsave(file.path(figdir, "Fig.2M.png"), plot=plot, width=3.5, height=3, dpi=300)

# Figure 2F: Pathway enrichment analysis of cell death pathways in β cell cluster of human islets exposed to mock or SARS-CoV-2 (MOI=1).
plot <- my.custom.bar.plot(infodir, suffix='DE.CoV2+.vs.Mock.beta_cells.mincells_10.wilcox')
ggsave(file.path(figdir, 'Fig.2F.png',sep='.'), width=6.5, height=5, dpi=300)

# Figure 6C: Dot plot analysis of the expression level of TNFSF12 in the macrophage cluster of human islets exposed to mock or SARS-CoV-2 virus (MOI=1).
# Load Seurat object with combined immune sub-cluster labels
panc.combined_labels <- readRDS(file.path(infodir,'panc.combined_labels.rds'))

# subset beta + macrophage cell clusters
tcells <- rownames(FetchData(panc.combined_labels, vars=c('ident','Name')) %>% 
  dplyr::filter(Name %in% c('control-9-mock','control-9-SARS-CoV-2','control-10-mock','control-10-SARS-CoV-2','control-11-mock','control-11-SARS-CoV-2')) %>% 
  dplyr::filter (ident %in% c(0,8)))

panc.cov2.bm <- subset(panc, cells=tcells)

tdata <- FetchData(panc.cov2.bm, vars=c('ident','Condition')) %>% 
  mutate(compare=paste0('C',ident,'-',Condition)) %>% dplyr::select('compare')

panc.cov2.bm %<>% AddMetaData(tdata)

plot <- DotPlot(subset(panc.cov2.bm, idents=c('C8-Mock','C8-COV2')), features=c('TNFSF12')) + 
  coord_flip() + theme(axis.text.x=element_text(angle=45, hjust=1))
ggsave(file.path(figdir, 'Fig.6C.png'), width=5, height=3.5, dpi=300)

# Figure 6B: Dot plot showed the differential signaling from macrophages to β cells in human islets exposed to mock or CVB4 virus 
cellchat.cvb4.merged <- readRDS(file.path(infodir, "cellchat.cvb4.merged.rds"))

png(file.path(figdir, 'Fig.6B.png'), width=4, height=3.5, units='in', res=300)
netVisual_bubble(cellchat.cvb4.merged, sources.use=c(7), targets.use = c(3),  comparison = c(1, 2), angle.x = 45, line.on=F)
dev.off()
