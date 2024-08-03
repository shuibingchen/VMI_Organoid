# my_functions.R
# customized functions for processing data and plotting
# 

library(Seurat, quietly=T, warn.conflicts=F)
library(dplyr, quietly=T, warn.conflicts=F)
library(ggplot2, quietly=T, warn.conflicts=F)
library(Matrix, quietly=T, warn.conflicts=F)
#library(scater, quietly=T, warn.conflicts=F)
library(CellChat, quietly=T, warn.conflicts=F)
library(reshape2, quietly=T, warn.conflicts=F)
library(pheatmap, quietly=T, warn.conflicts=F)
library(magrittr, quietly=T, warn.conflicts=F)
library(cowplot, quietly=T, warn.conflicts=F)
library(R.utils, quietly=T, warn.conflicts=F)
library(SingleCellExperiment, quietly=T, warn.conflicts=F)
library(tidyr, quietly=T, warn.conflicts=F)
library(tibble, quietly=T, warn.conflicts=F)
library(ggrepel, quietly=T, warn.conflicts=F)
library(ggplotify, quietly=T, warn.conflicts=F)

# read in raw counts data from a given sample
my.Read10X <- function(tcountdir, tprefix=NULL){
  # directory exists?
  if (! dir.exists(tcountdir)){
    print(paste("input raw counts folder does NOT exist:", tcountdir, sep=" "))
    return(NULL)
  }
  # file exists?
  tmat.file <- paste(tcountdir, "matrix.mtx.gz", sep="/")
  tfnames.file <- paste(tcountdir, "features.tsv.gz", sep="/")
  tbnames.file <- paste(tcountdir, "barcodes.tsv.gz", sep="/")
  for (tf in c(tmat.file, tfnames.file, tbnames.file)){
    if (! file.exists(tf)){
      print(paste("input file does NOT exist:", tf))
      return(NULL)
    }
  }
  # extract counts matrix
  cat(paste("Loading UMI counts table from", tcountdir, "..."))
  tmat <- readMM(paste(tcountdir, "matrix.mtx.gz", sep="/"))
  tfnames <- read.delim(paste(tcountdir, "features.tsv.gz", sep="/"), header=FALSE, stringsAsFactors=FALSE)
  tbnames <- read.delim(paste(tcountdir, "barcodes.tsv.gz", sep="/"), header=FALSE, stringsAsFactors=FALSE)
  # update column names (cell ids)
  if (is.null(tprefix)){
    colnames(tmat) = tbnames$V1
  }
  else{
    colnames(tmat) = paste(tprefix, tbnames$V1, sep="_")
  }
  #rownames(tmat) = tfnames$V1
  # replace rowname (Ensembl id) by gene symbol
  # in case gene symbol is not unique, append the _EnsemblID after it
  # missing gene symbol will be replaced by EnsemblID
  #tsymbols <- mapIds(EnsDb.Mmusculus.v79, keys=rownames(tmat), column="GENENAME", keytype="GENEID")
  ##tsymbols <- mapIds(EnsDb.Mmusculus.v79, keys=rownames(tmat), column="SYMBOL", keytype="GENEID")
  rownames(tmat) <- uniquifyFeatureNames(ID=tfnames$V1, names=tfnames$V2)
  cat(" done.","\n")
  return(tmat)
}

# merge raw read counts table collected from multiple samples, in a more efficient way
my.MergeMatrix.v2 <- function(tmats){
  cat("Merge raw UMI counts ")
  tfunc <- function(x,y){
    tres <- cbind(x,y[rownames(x),])
    cat(".")
    return(tres)
  }
  tmerged <- Reduce(f=tfunc, x=tmats)
  # fill na with 0
  tmerged[is.na(tmerged)] <- 0
  cat(" done.")
  return(tmerged)
}

# function to compare pre-defined two groups of cells
# select top DE genes and make volcano plot
my.DE.pair <- function(tobj, cdtA, cdtB, ntop, nCells, ttitle, no.mito=TRUE, mincells=10,
                       min.pct=0.1, logfc.threshold=0.1, test.use='wilcox',
                       cut.padj=0.1, cut.avglogfc=0,
                       cutFC=0.6, cutP=20, tfigdir, tinfodir, tsvg=FALSE){
  # check number of cells in each condition
  nA <- as.numeric(nCells %>% filter(ident==cdtA) %>% dplyr::select(n))
  nB <- as.numeric(nCells %>% filter(ident==cdtB) %>% dplyr::select(n))
  print(paste0('Compare ',cdtA,' (',nA,' cells)',' to ',cdtB, ' (',nB,' cells).'))
  # enough cells?
  if (nA < mincells | nB < mincells){
    print('Too few cells to perform DE.')
    return(NULL)
  }
  # run DE
  de.results <- FindMarkers(tobj, ident.1=cdtA, ident.2=cdtB, min.pct=min.pct, logfc.threshold=logfc.threshold, test.use=test.use)
  # write to file
  out.file.prefix <- paste('DE',cdtA,'vs',cdtB,paste0('mincells_',mincells),
                           test.use,paste0('min_pct_',min.pct),paste0('logfc_',logfc.threshold), sep='.')
  de.file <- file.path(tinfodir, paste(out.file.prefix, 'txt', sep='.'))
  write.table(de.results, de.file, quote=F, sep='\t', col.names=NA)
  # select top 10 up/down genes to label in the plot
  up.genes <- rownames(subset(de.results, p_val_adj < cut.padj & avg_log2FC > cut.avglogfc))
  down.genes <- rownames(subset(de.results, p_val_adj < cut.padj & avg_log2FC < -cut.avglogfc))
  top.up.genes <- head(up.genes, ntop)
  top.down.genes <- head(down.genes, ntop)
  # ignore mitochondiral genes?
  if (no.mito){
    top.up.genes <- head(grep(mito.pattern, up.genes, value=T, invert=T), ntop)
    top.down.genes <- head(grep(mito.pattern, down.genes, value=T, invert=T), ntop)
  }
  # make volcano plot
  g <- myVolcanoPlot(de.file, tx='p_val_adj', ty='avg_log2FC', tcutFC=cutFC, tcutP=cutP, tlabel.genes.up=top.up.genes, tlabel.genes.down=top.down.genes,
                     txlabel='LogFoldChange', tylabel='-LogP-value', tupper=300, talpha=0.8,
                     tcolor.up='firebrick3', tcolor.down='steelblue3', tcolor.other='gray60')
  g <- g + ggtitle(ttitle) + theme(plot.title=element_text(color='black',face='bold',size=20,hjust=0.5))
  ggsave(file.path(tfigdir,paste('Volcano',out.file.prefix,'png',sep='.')), plot=g, width=8, height=7, dpi=300)
  if (tsvg){
    ggsave(file.path(tfigdir,paste('Volcano',out.file.prefix,'svg',sep='.')), plot=g, width=8, height=7)
  }
  return(de.results)
}

# make volcano plot on a set of DE genes
myVolcanoPlot <- function(tdefile, tx='p_val', ty='avg_logFC', tcutFC=0.25, tcutP=20, tlabel.genes.up=NULL, tlabel.genes.down=NULL, txlabel='LogFoldChange', tylabel='-LogP-value', tupper=300, talpha=0.8, tcolor.up='firebrick3', tcolor.down='steelblue3', tcolor.other='gray60'){
  # read in DE results
  tde <- read.table(tdefile, sep='\t', header=T, check.names=F, stringsAsFactors=F, row.names=1)
  # arrange data
  tdataToPlot <- tde[,c(tx,ty)]
  tdataToPlot$gene <- rownames(tde)
  colnames(tdataToPlot) <- c('pval','logFC','gene')
  # log transform p-value
  tdataToPlot$logPval <- -log10(tdataToPlot$pval)
  # label and color
  if (is.null(tlabel.genes.up) | is.null(tlabel.genes.down)){
    tdataToPlot$color <- with(tdataToPlot, ifelse(logPval > tcutP & logFC > tcutFC, 'up', ifelse(logPval > tcutP & logFC < -tcutFC, 'down', 'other')))
    tdataToPlot$label <- with(tdataToPlot, ifelse(color %in% c('up','down'), gene, ''))
  } else{
    tdataToPlot$color <- with(tdataToPlot, ifelse(gene %in% tlabel.genes.up, 'up', ifelse(gene %in% tlabel.genes.down, 'down', 'other')))
    tdataToPlot$label <- with(tdataToPlot, ifelse(gene %in% c(tlabel.genes.up, tlabel.genes.down), gene, ''))
  }
  # any 0 p-values thus Inf logPval? modify to the upperlimit
  tdataToPlot$logPval[tdataToPlot$logPval > tupper] <- tupper
  # plot
  tg <- ggplot(tdataToPlot, aes(x=logFC, y=logPval, color=color, label=label))
  tg <- tg + geom_point(shape=19, size=2, alpha=talpha)
  tg <- tg + scale_color_manual(values=c('up'=tcolor.up,'down'=tcolor.down, 'other'=tcolor.other))
  tg <- tg + geom_text_repel()
  #tg <- tg + geom_vline(xintercept=tcutFC, linetype='dotted', color='gray75') + geom_vline(xintercept=-tcutFC, linetype='dotted', color='gray75')
  #tg <- tg + geom_hline(yintercept=tcutP, linetype="dotted", color="gray75")
  tg <- tg + xlab(txlabel) + ylab(tylabel)
  tg <- tg + theme_classic()
  tg <- tg + theme(legend.position='none')
  tg <- tg + theme(axis.title=element_text(size=18, color='black'), axis.text=element_text(size=16, color='black'))
  return(tg)
}

# pre-processing cells in a particular group
my.preprocessing.cellchat <- function(seurat.obj, sample.use, clusters.use, metadata, group.by='labels',
                                      CellChatDB.use=CellChatDB.use, population.size=TRUE, infodir=infodir, 
                                      suffix='', sources.use = c(3,7), targets.use = c(3,7), workers=4){
  # define cells used for running CellChat
  tcell.use <- rownames(subset(FetchData(seurat.obj, vars=c('ident', 'Name')), 
                               Name %in% sample.use & ident %in% clusters.use))
  
  # prepare input data for running CellChat
  tdata.input <- seurat.obj[['RNA']]@data[, tcell.use]
  tmeta <- metadata[tcell.use, , drop=F]
  
  # create a CellChat object
  cellchat <- createCellChat(object=tdata.input, meta=tmeta, group.by=group.by)
  print(cellchat)
  
  # set the used database in the object
  cellchat@DB <- CellChatDB.use
  
  # Preprocessing the expression data for cell-cell communication analysis
  # subset the expression data of signaling genes for saving computation cost
  print('subsetData')
  cellchat <- subsetData(cellchat) # This step is necessary even if using the whole database

  print('identifyOverExpressedGenes')
  cellchat <- identifyOverExpressedGenes(cellchat)
  print('identifyOverExpressedInteractions')
  cellchat <- identifyOverExpressedInteractions(cellchat)
  
  # project gene expression data onto PPI (Optional: when running it, USER should set `raw.use = FALSE` in the function `computeCommunProb()` in order to use the projected data)
  # cellchat <- projectData(cellchat, PPI.human)
  
  # Compute the communication probability and infer cellular communication network
  print('computeCommunProb')
  cellchat <- computeCommunProb(cellchat, population.size=population.size)
  # Filter out the cell-cell communication if there are only few number of cells in certain cell groups
  cellchat <- filterCommunication(cellchat, min.cells = 10)
  
  # Extract the inferred cellular communication network as a data frame
  df.net <- subsetCommunication(cellchat)
  write.table(df.net, file=file.path(infodir, paste('cell-cell.communication',suffix,'all.txt',sep='.')), quote=F, sep='\t', row.names=F)
  
  # interactions between beta and macrophage
  print(subsetCommunication(cellchat, sources.use = sources.use, targets.use = targets.use))
  
  # Infer the cell-cell communication at a signaling pathway level
  # CellChat computes the communication probability on signaling pathway level by summarizing the communication probabilities of all ligands-receptors interactions associated with each signaling pathway.
  # NB: The inferred intercellular communication network of each ligand-receptor pair and each signaling pathway is stored in the slot ‘net’ and ‘netP’, respectively.
  print('computeCommunProbPathway')
  cellchat <- computeCommunProbPathway(cellchat)
  
  # Calculate the aggregated cell-cell communication network
  print('aggregateNet')
  cellchat <- aggregateNet(cellchat)
  
  # Compute the network centrality scores
  print('netAnalysis_computeCentrality')
  future::plan("sequential")
  cellchat <- netAnalysis_computeCentrality(cellchat, slot.name = "netP")
  future::plan("multicore", workers = workers)
  
  # save CellChat object
  saveRDS(cellchat, file.path(infodir, paste('cellchat',suffix,'rds',sep='.')))
  
  return(cellchat)
}

