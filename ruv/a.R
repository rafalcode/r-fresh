#!/usr/bin/env Rscript
library(Cairo)
library(RUVSeq)
library(zebrafishRNASeq)
library(RColorBrewer)

# set colors early ...
colors <- brewer.pal(3, "Set2")
# introducing our data,
# 6 samples with genes counts.
# a typical dat matrix (here, in df form) where the rows are genes and there are 1000s of them
# an there are 6 columns, one for each sample.
data(zfGenes)
# get indices of the genes (rows) that have over 5 counts in at least two of the 6 samples.
filter <- apply(zfGenes, 1, function(x) length(x[x>5])>=2)
filtered <- zfGenes[filter,]

# get names of two known categories of genes.
genes <- rownames(filtered)[grep("^ENS", rownames(filtered))]
spikes <- rownames(filtered)[grep("^ERCC", rownames(filtered))]

# first phenotypical data, the 6 samples can be grouped, first three controls, second three treated's.
x <- as.factor(rep(c("Ctl", "Trt"), each=3))
# put'em together in a widelyused structure: expression set
# Expressionset belongs to Biobase package, but this function is from EDASeq.
set <- newSeqExpressionSet(as.matrix(filtered), phenoData = data.frame(x, row.names=colnames(filtered)))
set2 <- betweenLaneNormalization(set, which="upper")

CairoPNG("01_whiskplot_norming.png", 1600, 800)
par(mfrow=c(1,2))
plotRLE(set, outline=FALSE, ylim=c(-4, 4), col=colors[x])
plotRLE(set2, outline=FALSE, ylim=c(-4, 4), col=colors[x])
dev.off()


CairoPNG("02_pcaplot_norming.png", 1600, 800)
par(mfrow=c(1,2))
plotPCA(set, col=colors[x], cex=1.2)
plotPCA(set2, col=colors[x], cex=1.2)
dev.off()

diffs <- makeGroups(x)
# There are 3 key functions in RUVSeq packages: RUVg, RUVs, RUVr .. 
set3 <- RUVs(set2, genes, k=1, diffs)

CairoPNG("03_pcaplot_ruvsing.png", 1600, 800)
par(mfrow=c(1,2))
plotPCA(set2, col=colors[x], cex=1.2)
plotPCA(set3, col=colors[x], cex=1.2)
dev.off()

set4 <- RUVs(set2, genes, k=3, diffs)

CairoPNG("04_pcaplot_ruvsk3.png", 1600, 800)
par(mfrow=c(1,2))
plotPCA(set2, col=colors[x], cex=1.2)
plotPCA(set4, col=colors[x], cex=1.2)
dev.off()
