# my_functions.R
# customized functions for processing data and plotting
# 

library(Seurat, quietly=T, warn.conflicts=F)
library(Signac, quietly=T, warn.conflicts=F)
library(EnsDb.Hsapiens.v86)
library(tidyverse)
library(patchwork)
library(future)
library(magrittr)
library(RColorBrewer)
library(pheatmap)
library(GenomicRanges)

# ATAC blacklist regions
blacklist.regions <- blacklist_hg38_unified

# calculate cell QC metrics
my.compute.QC.metrics <- function(object){
    # compute nucleosome signal score per cell
    object %<>% NucleosomeSignal()
    # compute TSS enrichment score per cell
    object %<>% TSSEnrichment(fast = FALSE)
    # add fraction of reads in peaks
    object$pct_reads_in_peaks <- object$peak_region_fragments / object$passed_filters * 100
    ##object$blacklist_ratio <- object$blacklist_region_fragments / object$peak_region_fragments
    # compute blacklist ratio
    object$blacklist_fraction <- FractionCountsInRegion(object=object, library(Seurat, quietly=T, warn.conflicts=F)
                                                        assay='ATAC', 
                                                        regions=blacklist.regions)
    return(object)
}

# filter cells based on QC metrics
my.filter.cells <- function(object, nCount.min=3000, nCount.max=30000, pct.rip.min=20, 
                            bl.frac.max=0.05, nuc.sig.max=4, tss.min=3){
    return(subset(
        x = object,
        subset = nCount_ATAC > nCount.min & 
            nCount_ATAC < nCount.max &
            pct_reads_in_peaks > pct.rip.min &
            blacklist_fraction < bl.frac.max &
            nucleosome_signal < nuc.sig.max &
            TSS.enrichment > tss.min
    ))
}

# compute LSI
my.preprocess <- function(object, min.cutoff=20){
    # feature selection
    # include features with >20 total counts in the set of VariableFeatures
    object %<>% FindTopFeatures(min.cutoff = min.cutoff)
    # normalization
    object %<>% RunTFIDF()
    # dimension reduction
    object %<>% RunSVD()
    return(object)
}

myDimPlot3 <- function(tobj, treduct, tgroup_by, tgroup_order=NULL, thighlight=NULL, tsuffix, tcells=NULL, tcolor=NULL, tlabel=FALSE, tsplit_by=NULL, tsplit_order=NULL, txlim=NULL, tylim=NULL,
                      tncol=1, tptshape=19, tptsize=2, talpha=0.6, tltsize=18, tatlsize=20, tlbsize=2){
  # set coordinates variable name
  vars.reduct <- c("UMAP_1","UMAP_2")
  if (treduct == "tsne"){
    vars.reduct <- c("tSNE_1","tSNE_2")
  }
  # extract coordinates + group
  tdataToPlot <- FetchData(tobj, cells=tcells, vars=c(vars.reduct, tgroup_by))
  colnames(tdataToPlot) <- c("Dim_1","Dim_2","Group")
  # update group order if available
  if (!is.null(tgroup_order)){
    # add fake rows to make sure each group is considered
    tmiss.cates <- setdiff(tgroup_order, unique(tdataToPlot$Group))
    if (length(tmiss.cates) > 0){
      tdataToPlot <- rbind(tdataToPlot, data.frame(Dim_1=NA, Dim_2=NA, Group=tmiss.cates))
    }
    # reorder categories
    tdataToPlot$Group <- factor(tdataToPlot$Group, levels=tgroup_order)
  }
  # extract split
  if (! is.null(tsplit_by)){
    tsp <- FetchData(tobj, cells=tcells, vars=c(tsplit_by))
    tdataToPlot$Split <- tsp[rownames(tdataToPlot), c(tsplit_by)]
    # update split order if available
    if (! is.null(tsplit_order)){
      tdataToPlot$Split <- factor(tdataToPlot$Split, levels=tsplit_order)
    }
  }
  # reorder cells that needs to highlight (draw those cell points later)
  if (!is.null(thighlight)){
    tdataToPlot <- rbind(subset(tdataToPlot, ! Group %in% thighlight), subset(tdataToPlot, Group %in% thighlight))
  }
  # prepare group labeling
  tlabel.pos <- aggregate(cbind(Dim_1, Dim_2) ~ Group, data=tdataToPlot, FUN=median)
  colnames(tlabel.pos) <- c("Group","X","Y")
  # plot
  tg <- ggplot(tdataToPlot, aes(x=Dim_1, y=Dim_2, color=Group))
  tg <- tg + ggtitle(paste(toupper(treduct),"plots","by",tsuffix, sep=" "))
  # set range on coordinates
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
  if (! is.null(tsplit_by)){
    tg <- tg + facet_wrap(~Split, ncol=tncol)
  }
  if (tlabel == TRUE){
    tg <- tg + geom_text(data=tlabel.pos,aes(x=X, y=Y, label=Group), color="black", size=tlbsize)
  }
  return(tg)
}


