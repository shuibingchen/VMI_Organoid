# step_3.sample_integration.R
# --- process scRNAseq data from islets samples ---
# step 3: Integrate samples using MNN based corrections
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
library(sceasy)
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

# read in dissociation-related genes
disso.file <- file.path(sourcedata, "dissociation", "dissociation_related_genes.human.txt")
disso.genes <- as.vector(read.table(disso.file, header=F, check.names=F, sep="\t")$V1)

# sample info
sample.info <- data.frame(Name=c("control-9-mock","control-9-SARS-CoV-2","control-10-mock","control-10-SARS-CoV-2","control-10-CVB4","control-11-mock","control-11-SARS-CoV-2","control-12-mock","control-12-CVB4-1","control-12-CVB4-2"), 
                          Condition=c('Mock','COV2','Mock','COV2','CVB4','Mock','COV2','Mock','CVB4','CVB4'), 
                          Donor=c('control-9','control-9','control-10','control-10','control-10','control-11','control-11','control-12','control-12','control-12'),
                          Time=c('24h','24h','24h','24h','24h','24h','24h','24h','24h','48h'))
rownames(sample.info) <- c("D1S1","D1S2","D2S1","D2S2","D2S3","D3S1","D3S2","D4S1","D4S2","D4S3")

# load Seurat object
panc.initial <- readRDS(file.path(infodir, "panc.initial.rds"))

# --------------------------------------------- collect doublets --------------------------------------------- #
# collect and combine doublets info
stat <- data.frame(sample=character(), nTotal=integer(), nIsDoublet=integer(), nPreds=integer())
solo.res <- data.frame(solo.preds=character(), solo.isdoublet=character(), solo.score=numeric())
for (sid in unique(panc$orig.ident)){
  print(paste("processing",sid))
  # load cell ids
  cell.ids <- rownames(subset(FetchData(panc, vars=c('orig.ident')), orig.ident==sid))
  # read in solo results
  preds <- as.numeric(as.vector(read.table(file.path(solo.infodir, sid, 'preds.csv'), header=F, check.names=F, stringsAsFactors=F, sep=',')$V1))
  isdoublet <- as.numeric(as.vector(read.table(file.path(solo.infodir, sid, 'is_doublet.csv'), header=F, check.names=F, stringsAsFactors=F, sep=',')$V1))
  scores <- as.numeric(as.vector(read.table(file.path(solo.infodir, sid, 'logit_scores.csv'), header=F, check.names=F, stringsAsFactors=F, sep=',')$V1))
  # add to table
  data <- data.frame(solo.preds=ifelse(preds == 1, "doublet", "singlet"), solo.isdoublet=ifelse(isdoublet == 1, "doublet", "singlet"), solo.score=scores)
  rownames(data) <- cell.ids
  solo.res <- rbind(solo.res, data)
  # count doublets
  nPreds <- sum(data$solo.preds == "doublet")
  nIsDoublet <- sum(data$solo.isdoublet == "doublet")
  nTotal <- nrow(data)
  stat <- rbind(stat, data.frame(sample=sid, nTotal=c(nTotal), nIsDoublet=c(nIsDoublet), nPreds=c(nPreds)))
}

# save doublets result
write.table(solo.res, file=file.path(solo.infodir, 'solo.results.txt'), quote=F, sep='\t', row.names=T, col.names=NA)

# summarize doublets info
stat$pIsDoublet <- round(stat$nIsDoublet / stat$nTotal * 100, 2)
stat$pPreds <- round(stat$nPreds / stat$nTotal * 100, 2)
write.table(stat, file=file.path(solo.infodir, 'stats.solo.txt'), quote=F, sep='\t', row.names=F, col.names=T)

# --------------------------------------------- sample integration --------------------------------------------- #
# prepare raw UMI counts table from each donor
selected.donors <- c("D1S1","D1S2","D2S1","D2S2","D2S3","D3S1","D3S2","D4S1","D4S2","D4S3")
sample.list <- list()
for (donor in selected.donors){
  sample.list[[donor]] <- panc.initial[["RNA"]]@counts[, rownames(subset(panc.initial@meta.data, orig.ident == donor))]
}

# create SingleCellExperiment object
sce.list <- list()
for (donor in names(sample.list)){
  sce.list[[donor]] <- SingleCellExperiment(list(counts=as.matrix(sample.list[[donor]])))
}

# run a pre-clustering to avoid pooling together very different cells
# normalization will be performed for cells within each cluster
preclust.list <- lapply(sce.list, function(x) quickCluster(x=x, min.size=200, assay.type="counts", method="hclust", min.mean=0.1))

# normalize data by deconvolving size factors from cell pools
sce.list <- mapply(FUN=function(x,y) {computeSumFactors(x=x, min.mean=0.1, cluster=y)}, x=sce.list, y=preclust.list)

# compute normalized log-expression values
sce.list %<>% lapply(FUN=function(x) {normalize(object=x)})

# rescale among donors
rescaled.sce.list <- do.call(multiBatchNorm, sce.list)

# remove doublets
rescaled.sce.list <- lapply(rescaled.sce.list, FUN=function(x) { x[, intersect(colnames(x), rownames(subset(solo.res, solo.isdoublet=='singlet')))] })

# create a seurat object with raw UMI counts
panc <- CreateSeuratObject(counts=as(do.call(cbind, lapply(rescaled.sce.list, function(x) counts(x))), "dgCMatrix"), 
                           project=project, assay="RNA", min.cells=0, min.features=0,
                           names.field=1, names.delim="_", meta.data=NULL)

# Calculates the mitochondrial/ribosomal/viral genes per cell
panc[["percent.mito"]] <- PercentageFeatureSet(panc, pattern=mito.pattern)
panc[["percent.ribo"]] <- PercentageFeatureSet(panc, pattern=ribo.pattern)
panc[["percent.cov2"]] <- PercentageFeatureSet(panc, pattern=cov2.pattern)
panc[["percent.cvb4"]] <- PercentageFeatureSet(panc, pattern=cvb4.pattern)

# Add sample condition
tmeta <- data.frame(row.names=rownames(panc@meta.data))
for (tx in colnames(sample.info)){
  tdic <- as.vector(sample.info[,tx])
  names(tdic) <- rownames(sample.info)
  tmeta[,tx] <- as.vector(tdic[as.vector(panc@meta.data[,"orig.ident"])])
}
panc %<>% AddMetaData(metadata=tmeta)

# replace normalized data with the scran normalized data
panc[["RNA"]]@data <- as(do.call(cbind, lapply(rescaled.sce.list, function(x) logcounts(x))) * log(2), "dgCMatrix")

# Identification of highly variable features (feature selection)
panc %<>% FindVariableFeatures(selection.method="vst", nfeatures=3500)

# remove dissociation-related genes and ribosomal genes from variable gene list
# and select the remaining top 3000 genes for MNN-based correction
variable.genes <- setdiff(VariableFeatures(panc), c(disso.genes, grep(mito.pattern, rownames(panc), value=T), 
    grep(ribo.pattern, rownames(panc), value=T), grep(cov2.pattern, rownames(panc), value=T),
    grep(cvb4.pattern, rownames(panc), value=T)))
variable.genes <- head(variable.genes, 3000)

# perform MMN-based correction
original <- lapply(rescaled.sce.list, function(x) {logcounts(x)[variable.genes,]})
mnn.out <- do.call(fastMNN, c(original, list(k=20, d=50, auto.merge=TRUE)))
# set column names
colnames(reducedDim(mnn.out)) = paste0("MNN_", 1:ncol(reducedDim(mnn.out)))

# add MNN correction results to Seurat object
panc[["mnn"]] <- CreateDimReducObject(embeddings=reducedDim(mnn.out)[rownames(panc@meta.data),], key="MNN_", assay=DefaultAssay(panc))

# save Seurat object
saveRDS(panc, file=file.path(infodir, "panc.rds"))
