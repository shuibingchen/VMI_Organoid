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
library(clusterProfiler, quietly=T, warn.conflicts=FALSE)

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

# modify the Seurat DotPlot function to order samples
PercentAbove <- function(x, threshold) {
  return(length(x = x[x > threshold]) / length(x = x))
}

DotPlot.2 <- function (object, assay = NULL, features, cols = c("lightgrey", "blue"), col.min = -2.5, col.max = 2.5,
                       dot.min = 0, dot.scale = 6, group.by = NULL, split.by = NULL, scale.by = "radius",
                       scale.min = NA, scale.max = NA, order.ids = NULL){
  if (! is.null(assay)){
    assay <- DefaultAssay(object = object)
  }
  #assay <- assay %||% DefaultAssay(object = object)
  DefaultAssay(object = object) <- assay
  scale.func <- switch(EXPR = scale.by, size = scale_size,
                       radius = scale_radius, stop("'scale.by' must be either 'size' or 'radius'"))
  data.features <- FetchData(object = object, vars = features)
  data.features$id <- if (is.null(x = group.by)) {
    Idents(object = object)
  }
  else {
    object[[group.by, drop = TRUE]]
  }
  if (!is.factor(x = data.features$id)) {
    data.features$id <- factor(x = data.features$id)
  }
  id.levels <- levels(x = data.features$id)
  data.features$id <- as.vector(x = data.features$id)
  if (!is.null(x = split.by)) {
    splits <- object[[split.by, drop = TRUE]]
    if (length(x = unique(x = splits)) > length(x = cols)) {
      stop("Not enought colors for the number of groups")
    }
    cols <- cols[1:length(x = unique(x = splits))]
    names(x = cols) <- unique(x = splits)
    data.features$id <- paste(data.features$id, splits, sep = "_")
    unique.splits <- unique(x = splits)
    id.levels <- paste0(rep(x = id.levels, each = length(x = unique.splits)),
                        "_", rep(x = unique(x = splits), times = length(x = id.levels)))
  }
  data.plot <- lapply(X = unique(x = data.features$id), FUN = function(ident) {
    data.use <- data.features[data.features$id == ident,
                              1:(ncol(x = data.features) - 1), drop = FALSE]
    avg.exp <- apply(X = data.use, MARGIN = 2, FUN = function(x) {
      return(mean(x = expm1(x = x)))
    })
    pct.exp <- apply(X = data.use, MARGIN = 2, FUN = PercentAbove,
                     threshold = 0)
    return(list(avg.exp = avg.exp, pct.exp = pct.exp))
  })
  names(x = data.plot) <- unique(x = data.features$id)
  data.plot <- lapply(X = names(x = data.plot), FUN = function(x) {
    data.use <- as.data.frame(x = data.plot[[x]])
    data.use$features.plot <- rownames(x = data.use)
    data.use$id <- x
    return(data.use)
  })
  data.plot <- do.call(what = "rbind", args = data.plot)
  if (!is.null(x = id.levels)) {
    data.plot$id <- factor(x = data.plot$id, levels = id.levels)
  }
  avg.exp.scaled <- sapply(X = unique(x = data.plot$features.plot),
                           FUN = function(x) {
                             data.use <- data.plot[data.plot$features.plot ==
                                                     x, "avg.exp"]
                             data.use <- scale(x = data.use)
                             data.use <- MinMax(data = data.use, min = col.min,
                                                max = col.max)
                             return(data.use)
                           })
  avg.exp.scaled <- as.vector(x = t(x = avg.exp.scaled))
  if (!is.null(x = split.by)) {
    avg.exp.scaled <- as.numeric(x = cut(x = avg.exp.scaled,
                                         breaks = 20))
  }
  data.plot$avg.exp.scaled <- avg.exp.scaled
  data.plot$features.plot <- factor(x = data.plot$features.plot,
                                    levels = rev(x = features))
  data.plot$pct.exp[data.plot$pct.exp < dot.min] <- NA
  data.plot$pct.exp <- data.plot$pct.exp * 100
  if (!is.null(x = split.by)) {
    splits.use <- vapply(X = strsplit(x = as.character(x = data.plot$id),
                                      split = "_"), FUN = "[[", FUN.VALUE = character(length = 1L),
                         2)
    data.plot$colors <- mapply(FUN = function(color, value) {
      return(colorRampPalette(colors = c("grey", color))(20)[value])
    }, color = cols[splits.use], value = avg.exp.scaled)
  }
  color.by <- ifelse(test = is.null(x = split.by), yes = "avg.exp.scaled",
                     no = "colors")
  if (!is.na(x = scale.min)) {
    data.plot[data.plot$pct.exp < scale.min, "pct.exp"] <- scale.min
  }
  if (!is.na(x = scale.max)) {
    data.plot[data.plot$pct.exp > scale.max, "pct.exp"] <- scale.max
  }
  if (!is.null(order.ids)){
    data.plot$id <- factor(data.plot$id, levels=order.ids)
  }
  plot <- ggplot(data = data.plot, mapping = aes_string(x = "features.plot",
                                                        y = "id")) + geom_point(mapping = aes_string(size = "pct.exp",
                                                                                                     color = color.by)) + scale.func(range = c(0, dot.scale),
                                                                                                                                     limits = c(scale.min, scale.max)) + theme(axis.title.x = element_blank(),
                                                                                                                                                                               axis.title.y = element_blank()) + guides(size = guide_legend(title = "Percent Expressed")) +
    labs(x = "Features", y = ifelse(test = is.null(x = split.by),
                                    yes = "Identity", no = "Split Identity")) + theme_cowplot()
  if (!is.null(x = split.by)) {
    plot <- plot + scale_color_identity()
  }
  else if (length(x = cols) == 1) {
    plot <- plot + scale_color_distiller(palette = cols)
  }
  else {
    plot <- plot + scale_color_gradient(low = cols[1], high = cols[2])
  }
  if (is.null(x = split.by)) {
    plot <- plot + guides(color = guide_colorbar(title = "Average Expression"))
  }
  return(plot)
}

# function to compare pre-defined two groups of cells
# select top DE genes and make volcano plot
my.DE.pair <- function(tobj, cdtA, cdtB, ntop, nCells, ttitle, tinfodir, tfigdir, tsuffix=NULL, no.mito=TRUE, mincells=10,
                       min.pct=0.1, logfc.threshold=0.1, test.use='wilcox'){
  # check number of cells in each condition
  nA <- as.numeric(nCells %>% filter(ident==cdtA) %>% select(n))
  nB <- as.numeric(nCells %>% filter(ident==cdtB) %>% select(n))
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
  if (! is.null(tsuffix)){
    out.file.prefix <- paste('DE',cdtA,'vs',cdtB,tsuffix,paste0('mincells_',mincells),
                             test.use,paste0('min_pct_',min.pct),paste0('logfc_',logfc.threshold), sep='.')
  }
  de.file <- file.path(tinfodir, paste(out.file.prefix, 'txt', sep='.'))
  write.table(de.results, de.file, quote=F, sep='\t', col.names=NA)
  return(de.results)
}

my.run.GSEA.custom <- function(de.file, outdir, suffix, custom.pathways=custom.pathways, 
                        minGSSize=10, maxGSSize=500, 
                        pvalueCutoff=0.05, pAdjustMethod='BH', nterms=10, verbose=TRUE){
  # load DE results
  de.data <- read.table(de.file, header=T, check.names=F, stringsAsFactors=F, sep='\t', row.names=1)
  print(dim(de.data))
  
  # prepare input for GSEA (by Symbol)
  ## remove Covid/CVB4 genes
  clean.de.data <- de.data %>% rownames_to_column('Symbol') %>% 
    dplyr::filter(! grepl('^CoV2-|^polyprotein', Symbol)) %>% arrange(desc(avg_logFC))
  symbol.genes.gsea <- clean.de.data %>% pull(avg_logFC)
  names(symbol.genes.gsea) <- clean.de.data %>% pull('Symbol')
  
  # run GSEA on custom pathways (output all pathways)
  gse.custom <- GSEA(geneList=symbol.genes.gsea,
                     TERM2GENE=custom.pathways, 
                     minGSSize=minGSSize, 
                     maxGSSize=maxGSSize, 
                     pvalueCutoff=1, 
                     pAdjustMethod=pAdjustMethod, 
                     seed = TRUE, 
                     verbose=verbose)
  # save
  saveRDS(gse.custom, file.path(outdir, paste('gsea', paste('custom', suffix, sep='.'), 'rds', sep='.')))
}

# bar plot for GSEA
my.custom.bar.plot <- function(outdir, suffix){
  # load gsea result on custom pathways
  gsea.custom <- readRDS(file.path(outdir, paste('gsea', 'custom', suffix, 'rds', sep='.')))
  # calculating GeneRatio and Sign
  gsea.custom.forplot <- gsea.custom@result %>% mutate(GeneRatio=length(str_split(core_enrichment,'/')) / setSize,
                                             Sign=sign(enrichmentScore)) %>% arrange(GeneRatio)
  gsea.custom.forplot$ID <- factor(gsea.custom.forplot$ID, levels=as.vector(gsea.custom.forplot$ID))
  # plotting
  g <- ggplot(gsea.custom.forplot, aes(x=ID, y=GeneRatio, fill=enrichmentScore)) + geom_bar(stat="identity") 
  g <- g + scale_fill_distiller(palette='RdBu')
  g <- g + coord_flip() + theme_bw() + theme(axis.title.y = element_blank())
  print(g)
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
  #print(cellchat)
  
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
  future::plan("multiprocess", workers = workers)
  
  # save CellChat object
  saveRDS(cellchat, file.path(infodir, paste('cellchat',suffix,'rds',sep='.')))
  
  return(cellchat)
}
