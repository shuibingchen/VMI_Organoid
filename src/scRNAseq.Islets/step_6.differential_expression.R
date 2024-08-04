# step_6.differential_expression.R
# --- process scRNAseq data from islets samples ---
# step 6: perform differential expression analysis
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

# DE between SARS-N+ v.s. Mock

# beta cells in SARS-CoV-2 infected samples
beta.cov2.cells <- FetchData(panc, vars=c('Name','ident')) %>% rownames_to_column('cellID') %>% 
  dplyr::filter(Name %in% c('control-9-SARS-CoV-2','control-10-SARS-CoV-2','control-11-SARS-CoV-2') & ident %in% (0)) %>% pull('cellID')

# beta cells in SARS-CoV-2 Mock/Infected samples
beta.mock.cov2.cells <- FetchData(panc, vars=c('Name','ident')) %>% rownames_to_column('cellID') %>% 
  dplyr::filter(Name %in% c('control-9-mock','control-9-SARS-CoV-2','control-10-mock','control-10-SARS-CoV-2','control-11-mock','control-11-SARS-CoV-2') & ident %in% (0)) %>% pull('cellID')

# subset Seurat object
panc.beta.mock.cov2 <- subset(panc, cells=beta.mock.cov2.cells)

# separate SARS-N+ beta cells and SARS-N- beta cells 
cov2.genes <- grep(cov2.pattern, rownames(panc), value=T)

beta.mock.status <- FetchData(panc.beta.mock.cov2, vars='Condition') %>% rownames_to_column('cellID') %>% 
  dplyr::filter(Condition == 'Mock') %>% 
  mutate(Infection='Mock') %>% dplyr::select(c('cellID','Infection'))

beta.cov2.status <- data.frame(counts=colSums(panc.beta.mock.cov2[['RNA']]@counts[cov2.genes, beta.cov2.cells])) %>% rownames_to_column('cellID') %>% 
  mutate(Infection=ifelse(counts > 0, 'CoV2+', 'CoV2-')) %>% 
  dplyr::select(c('cellID', 'Infection'))

beta.mock.cov2.status <- rbind(beta.mock.status, beta.cov2.status) %>% column_to_rownames('cellID')

# add Mock/SARS-N+/- labels to Seurat
panc.beta.mock.cov2 %<>% AddMetaData(metadata=beta.mock.cov2.status)

# define comparison groups
Idents(panc.beta.mock.cov2) <- 'Infection'

# get number of cells per group
nCells.cdt <- FetchData(panc.beta.mock.cov2, vars=('Infection')) %>% dplyr::count(Infection) %>% dplyr::rename('ident'='Infection')

# DE: SARS-N+ v.s. Mock
tde <- my.DE.pair(tobj=panc.beta.mock.cov2, cdtA='CoV2+', cdtB='Mock', ntop=20, nCells=nCells.cdt, ttitle='beta.CoV2+.vs.Mock',
                  tinfodir=infodir, tfigdir=figdir, tsuffix='beta_cells', no.mito=TRUE, mincells=10, 
                  min.pct=0, logfc.threshold=0, test.use='wilcox')

# set back ident
Idents(panc) <- 'final.clust'

# perform GSEA

# DE result file
de.file <- file.path(infodir, 'DE.CoV2+.vs.Mock.beta_cells.mincells_10.wilcox.min_pct_0.logfc_0.txt')

my.run.GSEA.custom(de.file=de.file, outdir=infodir, 
	suffix='DE.CoV2+.vs.Mock.beta_cells.mincells_10.wilcox', 
	custom.pathways=custom.pathways, 
	minGSSize=10, maxGSSize=1000, 
	pvalueCutoff=0.05, pAdjustMethod='BH', nterms=10, verbose=FALSE)

my.custom.bar.plot(outdir, suffix='DE.CoV2+.vs.Mock.beta_cells.mincells_10.wilcox')

