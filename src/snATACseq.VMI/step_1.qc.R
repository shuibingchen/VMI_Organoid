# step_1.qc.R
# --- process scRNAseq data from VMI samples ---
# step 1: Load ATAC peak data and perform cell QC
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

# ATAC blacklist regions
blacklist.regions <- blacklist_hg38_unified

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

# extract gene annotations from EnsDb
annotations <- GetGRangesFromEnsDb(ensdb = EnsDb.Hsapiens.v86)

# change to USCS style since the chromsomes start with 'chr'
seqlevels(annotations) <- ifelse(seqlevels(annotations) == 'MT', 'chrM', paste0('chr', seqlevels(annotations)))
genome(annotations) <- "hg38"

# Creating a common peak set
gr.list <- list()

for (k in 1:nrow(sample.info)){
    pid <- sample.info$Name[k]
    sid <- sample.info$SeqName[k]
    print(sid)
    peaks <- read.table(file=file.path(sourcedir, sid, 'peaks.bed'), 
                        col.names = c("chr", "start", "end"))
    gr.list[[pid]] <- makeGRangesFromDataFrame(peaks)
}

# Create a unified set of peaks to quantify in each dataset
combined.peaks <- Signac::reduce(c(gr.list[[1]], gr.list[[2]], gr.list[[3]]))

# Filter out bad peaks based on length
peakwidths <- width(combined.peaks)
combined.peaks <- combined.peaks[peakwidths  < 10000 & peakwidths > 20]

# Create Seurat objects
seurat.obj.list <- list()

for (k in 1:nrow(sample.info)){
    pid <- sample.info$Name[k]
    sid <- sample.info$SeqName[k]
    # load metadata
    metadata <- read.table(file=file.path(sourcedir, sid, 'singlecell.csv'), header=T, 
                           stringsAsFactors=F, sep=',', row.names=1)
    # remove the first row
    metadata <- metadata[-1, ]
    # fetch valid cells
    metadata <- metadata[metadata$is__cell_barcode == 1, ]
    cells <- rownames(metadata)
    # create fragment objects
    frags <- CreateFragmentObject(path=file.path(sourcedir, sid, 'fragments.tsv.gz'), 
                                             cells=cells)
    # quantify peaks
    # create a matrix of peaks x cell
    counts <- FeatureMatrix(fragments=frags, features=combined.peaks, cells=cells)
    # create object
    chrom.assay <- CreateChromatinAssay(counts=counts, fragments=frags)
    seurat.obj.list[[pid]] <- CreateSeuratObject(counts=chrom.assay, assay = "ATAC", meta.data=metadata)
}

# add the gene information to the object
for (pid in sample.info$Name){
    Annotation(seurat.obj.list[[pid]]) <- annotations
}

# calculate cell QC metrics
seurat.obj.list <- lapply(seurat.obj.list, FUN=my.compute.QC.metrics)

# filter cells based on QC metrics
seurat.obj.list <- lapply(seurat.obj.list, FUN=my.filter.cells, 
                          nCount.min=3000, nCount.max=30000, pct.rip.min=20, 
                          bl.frac.max=0.05, nuc.sig.max=4, tss.min=3)

# save seurat object
saveRDS(seurat.obj.list, file.path(infodir, 'seurat.obj.list.rds'))
