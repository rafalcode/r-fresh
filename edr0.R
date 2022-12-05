#!/usr/bin/env Rscript
# edgeR simulations of count matrice use rnbinom()
# actually edgeR is not much into simulations
# or just views them as simplistic.
# so they are only mentione dunder the plotMDS function.
library(Cairo)
library(edgeR)

# Simulate DGE data for 1000 genes and 6 samples.
# Samples are in two groups
# First 200 genes are differentially expressed in secon
ngenes <- 1000
nlib <- 6 # sample library
counts <- matrix(rnbinom(ngenes*nlib, size=1/10, mu=20), ngenes, nlib)
rownames(counts) <- sprintf("Gene%04i", 1:ngenes)
colnames(counts) <- sprintf("S%02i", 1:nlib)

group <- gl(2,3,labels=c("Grp1","Grp2"))
counts[1:200,group=="Grp2"] <- counts[1:200,group=="Grp2"] +10 # just add ten Grp2
y <- DGEList(counts,group=group)
y <- calcNormFactors(y)
col <- as.numeric(group) # factor to numeric, not usually done.

# without labels, indexes of samples are plotted.
CairoPNG("edrmds0.png", 800, 800)
mds <- plotMDS(y, top=200, col=col)
# mds <- plotMDS(y, col=col)
dev.off()

#with labs
CairoPNG("edrmds1.png", 800, 800)
# if you saved the returned datastruct from above, you can do an identical graph thisaway
plotMDS(mds, col=col, labels=group)
# plotMDS(y, col=col, labels=group)
dev.off()
