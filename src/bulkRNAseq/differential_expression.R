# differential_expression.R
# perform differential analysis on bulk RNA-seq samples using DESeq2
# Author: Tuo Zhang
# Date: 8/1/2024
# 

library(DESeq2)
library(ggplot2)
library(vsn)
library(pheatmap)
library(RColorBrewer)
library(gplots)
library(BiocParallel)
library(apeglm)

workdir <- "."
sourcedir <- file.path(workdir, "data")
figdir <- file.path(workdir, "figure")
numdir <- file.path(workdir, "info")

# false discovery rate cutoff
alpha <- 0.1

# use parallelization
register(MulticoreParam(4))

# read in raw counts file
counts.file <- file.path(sourcedir, "bulk", "raw_counts.txt", sep="/")
countData <- read.table(file=countsfile, header=TRUE, check.names=FALSE)
rownames(countData) <- countData[,1]
countData <- countData[,-1]

#M1-1    M1-2    M1-3    M0-1    M0-2    M0-3

# create experimental label
colData <- data.frame(condition=factor(c(rep("M1",3),rep("M0",3))))
rownames(colData) <- colnames(countData)

# construct a DESeqDataSet
dds <- DESeqDataSetFromMatrix(countData=countData, colData=colData, design=~condition)
dds

# update factor levels
dds$condition <- factor(dds$condition, levels=c("M0","M1"))

# run DESeq
dds <- DESeq(dds, parallel=TRUE, BPPARAM=MulticoreParam(4))

# collect DE genes
res <- results(dds, contrast=c("condition", 'M0', 'M1'), alpha=alpha, parallel=TRUE, BPPARAM=MulticoreParam(4))

# shrink log2 fold change values
resLFC <- lfcShrink(dds, res=res, type="ashr", parallel=TRUE, BPPARAM=MulticoreParam(4))
			
# reorder results by adjusted-pvalue
resLFCOrdered <- resLFC[order(resLFC$padj),]

# make MA-plot
png(file=paste(figdir, paste("MA-plot", allconditions[i], "vs", allconditions[j], "LFC", "png", sep="."), sep="/"), width = 1920, height = 1920, units = "px", pointsize = 6, res=600)
plotMA(resLFC, main=paste(allconditions[i], "vs", allconditions[j], sep=" "))
abline(h=c(-1,1), col="dodgerblue", lwd=2)
dev.off()

# exporting results to file
write.table(as.data.frame(resLFCOrdered), file=file.path(numdir, "DE.M1.vs.M0.LFC.tsv"), quote=FALSE, sep='\t', col.names=NA)
		
# data visualization
rld <- rlog(dds, blind=F)

# heatmap of the sample-to-sample distances
sampleDists <- dist(t(assay(rld)))
sampleDistMatrix <- as.matrix(sampleDists)
hmcol <- colorRampPalette(brewer.pal(9, "GnBu"))(255)
png(file=paste(figdir, "heatmap.sample-to-sample.rld.png", sep="/"), width = 3600, height = 3200, units = "px", pointsize = 6, res=600)
pheatmap(sampleDistMatrix, clustering_distance_rows=sampleDists, clustering_distance_cols=sampleDists, col=hmcol, fontsize_row=6, fontsize_col=6)
dev.off()

# updated PCA plot with sample labels
mypcadata <- plotPCA(rld, intgroup=c("condition"), returnData=TRUE)
percentVar <- round(100*attr(mypcadata, "percentVar"))
g <- ggplot(mypcadata, aes(PC1, PC2, color=condition, label=name))
g <- g + geom_point(size=3, shape=19) + xlab(paste0("PC1: ", percentVar[1], "% variance")) + ylab(paste0("PC2: ", percentVar[2], "% variance"))
g <- g + theme_classic()
ggsave(file=paste(figdir, "PCA.rld.v2.png", sep="/"), width = 8, height = 6, type = "cairo", dpi = 600)

g <- ggplot(mypcadata, aes(PC1, PC2, color=condition, label=name))
g <- g + geom_text() + xlab(paste0("PC1: ", percentVar[1], "% variance")) + ylab(paste0("PC2: ", percentVar[2], "% variance"))
g <- g + theme_classic()
ggsave(file=paste(figdir, "PCA.rld.v2.labeled.png", sep="/"), width = 8, height = 6, type = "cairo", dpi = 600)

# regularized log transformation
write.table(as.data.frame(assay(rld)), file=paste(numdir, "details_rld.tsv", sep="/"), quote=FALSE, sep='\t', col.names=NA)

# save R object
saveRDS(dds, file=paste(numdir, "dds.rds", sep="/"))
saveRDS(rld, file=paste(numdir, "rld.rds", sep="/"))

sessionInfo()
