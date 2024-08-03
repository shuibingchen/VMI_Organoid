# my_functions.R
# customized functions for processing data and plotting
# 

library(Seurat, quietly=T, warn.conflicts=F)
library(dplyr, quietly=T, warn.conflicts=F)
library(ggplot2, quietly=T, warn.conflicts=F)
library(Matrix, quietly=T, warn.conflicts=F)
library(scater, quietly=T, warn.conflicts=F)
library(reshape2, quietly=T, warn.conflicts=F)
library(pheatmap, quietly=T, warn.conflicts=F)
library(magrittr, quietly=T, warn.conflicts=F)
library(cowplot, quietly=T, warn.conflicts=F)
library(R.utils, quietly=T, warn.conflicts=F)

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

# plot cells colored by an annotation (DimPlot)
myDimPlot <- function(tobj, treduct, tcate, torder=NULL, tsuffix, tcells=NULL, tcolor=NULL, tlabel=FALSE, tsplit=FALSE, txlim=NULL, tylim=NULL,
                      tptshape=19, tptsize=2, talpha=0.6, tltsize=18, tatlsize=20){
  tdataToPlot <- data.frame()
  tlabel.pos <- data.frame()
  tg <- ggplot()
  if (treduct == "tsne"){
    if (is.null(tcells)){
      tdataToPlot <- FetchData(tobj, vars=c("tSNE_1","tSNE_2",tcate))
    } else {
      tdataToPlot <- FetchData(tobj, vars=c("tSNE_1","tSNE_2",tcate), cells=tcells)
    }
    colnames(tdataToPlot) <- c("tSNE_1","tSNE_2","Category")
    if (!is.null(torder)){
      # add fake rows to make sure each category is considered
      tdataToPlot <- rbind(tdataToPlot, data.frame(tSNE_1=NA, tSNE_2=NA, Category=torder))
      # reorder categories
      tdataToPlot$Category <- factor(tdataToPlot$Category, levels=torder)
    }
    tg <- ggplot(tdataToPlot, aes(x=tSNE_1, y=tSNE_2, color=Category))
    tg <- tg + ggtitle(paste("tSNE","plots","by",tsuffix, sep=" "))
    tlabel.pos <- aggregate(cbind(tSNE_1, tSNE_2) ~ Category, data=tdataToPlot, FUN=median)
    colnames(tlabel.pos) <- c("Category","X","Y")
  } else if (treduct == "umap") {
    if (is.null(tcells)){
      tdataToPlot <- FetchData(tobj, vars=c("UMAP_1","UMAP_2",tcate))
    } else {
      tdataToPlot <- FetchData(tobj, vars=c("UMAP_1","UMAP_2",tcate), cells=tcells)
    }
    colnames(tdataToPlot) <- c("UMAP_1","UMAP_2","Category")
    if (!is.null(torder)){
      # add fake rows to make sure each category is considered
      tdataToPlot <- rbind(tdataToPlot, data.frame(UMAP_1=NA, UMAP_2=NA, Category=torder))
      # reorder categories
      tdataToPlot$Category <- factor(tdataToPlot$Category, levels=torder)
    }
    tg <- ggplot(tdataToPlot, aes(x=UMAP_1, y=UMAP_2, color=Category))
    tg <- tg + ggtitle(paste("UMAP","plots","by",tsuffix, sep=" "))
    tlabel.pos <- aggregate(cbind(UMAP_1, UMAP_2) ~ Category, data=tdataToPlot, FUN=median)
    colnames(tlabel.pos) <- c("Category","X","Y")
  }
  if (! is.null(txlim)){
    tg <- tg + coord_cartesian(xlim=txlim, ylim=tylim)
  }
  tg <- tg + theme_bw() + theme(panel.grid.major=element_blank(), panel.grid.minor=element_blank(), legend.justification=c(0,0), legend.title=element_blank())
  tg <- tg + theme(legend.key=element_blank()) + theme(legend.text=element_text(size=tltsize))
  tg <- tg + theme(axis.text=element_blank(), axis.title=element_text(size=tatlsize,face="bold"))
  tg <- tg + theme(axis.ticks=element_blank())
  tg <- tg + theme(plot.title=element_text(hjust=0.5))
  tg <- tg + geom_point(shape=tptshape, size=tptsize, alpha=talpha)
  if (! is.null(tcolor)){
    tg <- tg + scale_color_manual(values=tcolor)
  }
  if (treduct == "tsne"){
    tg <- tg + xlab("tSNE 1") + ylab("tSNE 2")
  } else {
    tg <- tg + xlab("UMAP 1") + ylab("UMAP 2")
  }
  if (tsplit == TRUE){
    tg <- tg + facet_wrap(~Category)
  }
  if (tlabel == TRUE){
    tg <- tg + geom_text(data=tlabel.pos,aes(x=X, y=Y, label=Category), color="black")
  }
  return(tg)
}

# plot cells colored by an annotation (DimPlot) (line and fill)
myDimPlot2 <- function(tobj, treduct, tcate, torder=NULL, tsuffix, tcells=NULL, tcolor=NULL, tedgecolor='gray50', tlabel=FALSE, tsplit=FALSE, txlim=NULL, tylim=NULL,
                      tptshape=21, tptsize=2, tstroke=0.3, talpha=0.6, tltsize=18, tatlsize=20){
  tdataToPlot <- data.frame()
  tlabel.pos <- data.frame()
  tg <- ggplot()
  if (treduct == "tsne"){
    if (is.null(tcells)){
      tdataToPlot <- FetchData(tobj, vars=c("tSNE_1","tSNE_2",tcate))
    } else {
      tdataToPlot <- FetchData(tobj, vars=c("tSNE_1","tSNE_2",tcate), cells=tcells)
    }
    colnames(tdataToPlot) <- c("tSNE_1","tSNE_2","Category")
    if (!is.null(torder)){
      # add fake rows to make sure each category is considered
      tdataToPlot <- rbind(tdataToPlot, data.frame(tSNE_1=NA, tSNE_2=NA, Category=torder))
      # reorder categories
      tdataToPlot$Category <- factor(tdataToPlot$Category, levels=torder)
    }
    tg <- ggplot(tdataToPlot, aes(x=tSNE_1, y=tSNE_2, fill=Category))
    tg <- tg + ggtitle(paste("tSNE","plots","by",tsuffix, sep=" "))
    tlabel.pos <- aggregate(cbind(tSNE_1, tSNE_2) ~ Category, data=tdataToPlot, FUN=median)
    colnames(tlabel.pos) <- c("Category","X","Y")
  } else if (treduct == "umap") {
    if (is.null(tcells)){
      tdataToPlot <- FetchData(tobj, vars=c("UMAP_1","UMAP_2",tcate))
    } else {
      tdataToPlot <- FetchData(tobj, vars=c("UMAP_1","UMAP_2",tcate), cells=tcells)
    }
    colnames(tdataToPlot) <- c("UMAP_1","UMAP_2","Category")
    if (!is.null(torder)){
      # add fake rows to make sure each category is considered
      tdataToPlot <- rbind(tdataToPlot, data.frame(UMAP_1=NA, UMAP_2=NA, Category=torder))
      # reorder categories
      tdataToPlot$Category <- factor(tdataToPlot$Category, levels=torder)
    }
    tg <- ggplot(tdataToPlot, aes(x=UMAP_1, y=UMAP_2, fill=Category))
    tg <- tg + ggtitle(paste("UMAP","plots","by",tsuffix, sep=" "))
    tlabel.pos <- aggregate(cbind(UMAP_1, UMAP_2) ~ Category, data=tdataToPlot, FUN=median)
    colnames(tlabel.pos) <- c("Category","X","Y")
  }
  if (! is.null(txlim)){
    tg <- tg + coord_cartesian(xlim=txlim, ylim=tylim)
  }
  tg <- tg + theme_bw() + theme(panel.grid.major=element_blank(), panel.grid.minor=element_blank(), legend.justification=c(0,0), legend.title=element_blank())
  tg <- tg + theme(legend.key=element_blank()) + theme(legend.text=element_text(size=tltsize))
  tg <- tg + theme(axis.text=element_blank(), axis.title=element_text(size=tatlsize,face="bold"))
  tg <- tg + theme(axis.ticks=element_blank())
  tg <- tg + theme(plot.title=element_text(hjust=0.5))
  tg <- tg + geom_point(shape=tptshape, size=tptsize, alpha=talpha, color=tedgecolor, stroke=tstroke)
  if (! is.null(tcolor)){
    tg <- tg + scale_fill_manual(values=tcolor)
  }
  if (treduct == "tsne"){
    tg <- tg + xlab("tSNE 1") + ylab("tSNE 2")
  } else {
    tg <- tg + xlab("UMAP 1") + ylab("UMAP 2")
  }
  if (tsplit == TRUE){
    tg <- tg + facet_wrap(~Category)
  }
  if (tlabel == TRUE){
    tg <- tg + geom_text(data=tlabel.pos,aes(x=X, y=Y, label=Category), color="black")
  }
  return(tg)
}

# highlight expression of a set of genes (FeaturePlot)
MyFeaturePlot <- function(tobj, tgenes, tcells=NULL, tassay="RNA", treduction.name="umap", txlim=NULL, tylim=NULL, tbreaks=NULL, tlimits=NULL, tlowcolor='gray80', thighcolor='red2', tncol=2, tlegend=NULL, tptsize=2, talpha=0.7, twidth=15, theight=12.5, tunits="in", tres=300){
  # genes valid?
  tgenes.valid <- intersect(tgenes, rownames(tobj))
  if (is.null(tgenes.valid)){
    cat("No valid genes found, do nothing!")
    return(NULL)
  }
  # assay valid?
  if (! tassay %in% names(tobj)){
    cat(paste("Not a valid assay:",tassay,sep=" "))
    return(NULL)
  }
  # extract gene expression
  texp <- as.matrix(tobj[[tassay]]@data[tgenes.valid, ,drop=F])
  # get coordinates
  tvars <- c("UMAP_1","UMAP_2")
  if (treduction.name == "tsne"){
    tvars <- c("tSNE_1","tSNE_2")
  }
  tdata <- FetchData(object=tobj, vars=tvars)
  colnames(tdata) <- c("X","Y")
  # plot
  tplots <- list()
  tk <- 1
  for (tgene in tgenes){
    # merge data for plotting
    tdata.merged <- merge(tdata, t(texp[tgene,,drop=F]), by=0, all=T)
    rownames(tdata.merged) <- tdata.merged$Row.names
    tdata.merged <- tdata.merged[,-1]
    colnames(tdata.merged) <- c("X","Y","Expression")
    # subset cells?
    if (! is.null(tcells)){
      tdata.merged <- tdata.merged[tcells,,drop=F]
    }
    # reorder cells by expression of the given gene
    tdata.merged <- tdata.merged[with(tdata.merged, order(Expression)),]
    # plot
    if (max(tdata.merged$Expression) > 0){ # expressed in at least one cell
      # plot (rename x and y axis)
      tg <- ggplot(tdata.merged, aes(x=X, y=Y, color=Expression))
      tg <- tg + geom_point(shape=19, size=tptsize, alpha=talpha)
      if (! is.null(tbreaks)){
        tg <- tg + scale_color_gradient(low=tlowcolor, high=thighcolor, breaks=tbreaks, limits=tlimits)
      } else {
        tg <- tg + scale_color_gradient(low=tlowcolor, high=thighcolor)
      }
      if(! is.null(txlim)){
        tg <- tg + coord_cartesian(xlim=txlim, ylim=tylim)
      }
      tg <- tg + ggtitle(tgene)
      if (treduction.name == "tsne"){
        tg <- tg + xlab("tSNE 1") + ylab("tSNE 2")
      } else {
        tg <- tg + xlab("UMAP 1") + ylab("UMAP 2")
      }
      tg <- tg + theme_bw()
      tg <- tg + theme(plot.title=element_text(hjust=0.5, size=18, face="bold"))
      tg <- tg + theme(axis.text=element_blank(), axis.ticks=element_blank(), axis.title=element_text(size=16,face="bold"))
      # add to list
      tplots[[tk]] <- tg
      tk <- tk + 1
    } else { # no expressions at all
      tg <- ggplot(tdata.merged, aes(x=X, y=Y))
      tg <- tg + geom_point(color="gray80", shape=19, size=tptsize, alpha=talpha)
      if(! is.null(txlim)){
        tg <- tg + coord_cartesian(xlim=txlim, ylim=tylim)
      }
      tg <- tg + ggtitle(tgene)
      if (treduction.name == "tsne"){
        tg <- tg + xlab("tSNE 1") + ylab("tSNE 2")
      } else {
        tg <- tg + xlab("UMAP 1") + ylab("UMAP 2")
      }
      tg <- tg + theme_bw()
      tg <- tg + theme(plot.title=element_text(hjust=0.5, size=18, face="bold"))
      tg <- tg + theme(axis.text=element_blank(), axis.ticks=element_blank(), axis.title=element_text(size=16,face="bold"))
      # add to list
      tplots[[tk]] <- tg
      tk <- tk + 1
    }
  }
  # combine plots with Seurat::CombinePlots
  tcombined <- CombinePlots(tplots, ncol=tncol, legend=tlegend)
  return(tcombined)
  # ggsave(paste(toutdir, paste("gene","exp",tsuffix,toupper(treduction.name),"png",sep="."),sep="/"), height=theight, width=twidth, units=tunits, dpi=tres)
  # ggsave(paste(toutdir, paste("gene","exp",tsuffix,toupper(treduction.name),"svg",sep="."),sep="/"), height=theight, width=twidth, units=tunits)
}


