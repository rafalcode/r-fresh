#!/usr/bin/env Rscript
# this script does what? After frsutration with my customised pcaplot ... I go back to drawing boaard
# this is taken from the DESeq refman.
# note: he does refer to the vignette which also has an example, using vst()

# Beware the use of "transforamtion" here. It's a sort of fondling of the data
# which helps alot with visual tools such as plots.

library(DESeq2)
library(ggplot2)
library(Cairo)
library(ggrepel)

plotPCAorig <- function(object, intgroup="condition", ntop=500, returnData=FALSE)
{
  # calculate the variance for each gene
  rv <- rowVars(assay(object))

  # select the ntop genes by variance
  select <- order(rv, decreasing=TRUE)[seq_len(min(ntop, length(rv)))]

  # perform a PCA on the data in assay(x) for the selected genes
  pca <- prcomp(t(assay(object)[select,]))

  # the contribution to the total variance for each component
  percentVar <- pca$sdev^2 / sum( pca$sdev^2 )

  if (!all(intgroup %in% names(colData(object)))) {
    stop("the argument 'intgroup' should specify columns of colData(dds)")
  }

  intgroup.df <- as.data.frame(colData(object)[, intgroup, drop=FALSE])
  
  # add the intgroup factors together to create a new grouping factor
  group <- if (length(intgroup) > 1) {
    factor(apply( intgroup.df, 1, paste, collapse=":"))
  } else {
    colData(object)[[intgroup]]
  }

  # assembly the data for the plot
  d <- data.frame(PC1=pca$x[,1], PC2=pca$x[,2], group=group, intgroup.df, name=colnames(object))

  if (returnData) {
    attr(d, "percentVar") <- percentVar[1:2]
    return(d)
  }
  
  ggplot(data=d, aes_string(x="PC1", y="PC2", color="group")) + geom_point(size=3) + 
    xlab(paste0("PC1: ",round(percentVar[1] * 100),"% variance")) +
      ylab(paste0("PC2: ",round(percentVar[2] * 100),"% variance")) +
        coord_fixed()
}

plotPCArf <- function(object, intgroup="condition", ntop=500, returnData=FALSE)
{
  # calculate the variance for each gene
  rv <- rowVars(assay(object))

  # select the ntop genes by variance
  select <- order(rv, decreasing=TRUE)[seq_len(min(ntop, length(rv)))]

  # perform a PCA on the data in assay(x) for the selected genes
  pca <- prcomp(t(assay(object)[select,]))

  # the contribution to the total variance for each component
  percentVar <- pca$sdev^2 / sum( pca$sdev^2 )

  if (!all(intgroup %in% names(colData(object)))) {
    stop("the argument 'intgroup' should specify columns of colData(dds)")
  }

  intgroup.df <- as.data.frame(colData(object)[, intgroup, drop=FALSE])
  
  # add the intgroup factors together to create a new grouping factor
  group <- if (length(intgroup) > 1) {
    factor(apply( intgroup.df, 1, paste, collapse=":"))
  } else {
    colData(object)[[intgroup]]
  }

  # assembly the data for the plot
  d <- data.frame(PC1=pca$x[,1], PC2=pca$x[,2], group=group, intgroup.df, name=colnames(object))

  if (returnData) {
    attr(d, "percentVar") <- percentVar[1:2]
    return(d)
  }
  
  ggplot(data=d, aes(x=!!sym("PC1"), y=!!sym("PC2"), color=!!sym("group"), label=!!sym("name"))) + geom_point(size=3) + 
  # ggplot(data=d, aes_string(x="PC1", y="PC2", color="group", label="name")) + geom_point(size=3) + 
    xlab(paste0("PC1: ",round(percentVar[1] * 100),"% variance")) +
      ylab(paste0("PC2: ",round(percentVar[2] * 100),"% variance")) +
       geom_text_repel(size=4, position = position_nudge_repel(x = 0.1, y = 0.1)) + coord_fixed()
} # end of funct

dds <- makeExampleDESeqDataSet(betaSD=1)
# by default (and you should know by now)
# you get 1k genes in 12 samples, 6 each per the 2 conditions.

# colnames(dds) <- gsub("[^0-9]+(0-9+)", "S\\1", colnames(dds))
# don;t know why that didn't work ... the following will:
colnames(dds) <- gsub("[^\\d]+(\\d+)$", "S\\1", colnames(dds), perl=T)

# using rlog transformed data:
# transform either with rlog or vst
rld <- rlog(dds)

# Cairo image template
CairoPNG("despca01.png", 800, 800)
# the way I called it:
# pp <- plotPCArf(vsdds, intgroup="Tissue") # my version of the function.
plotPCArf(rld)
# plotPCAorig(rld)
dev.off()

# He introdcues a second part for some reason:
# also possible to perform custom transformation:
# this time, there is no explicit transformation
# (no intermediate transofrmation variable created).
dds2 <- estimateSizeFactors(dds) # obviously necessary.

# shifted log of normalized counts , the shift may even jsut be the addaiton of 1.
se <- SummarizedExperiment(log2(counts(dds2, normalized=TRUE) + 1), colData=colData(dds2))

# the call to DESeqTransform() is needed to
# trigger our plotPCA method. That's another way of putting it (RF)
CairoPNG("despca02.png", 800, 800)
plotPCA( DESeqTransform(se) )
dev.off()

