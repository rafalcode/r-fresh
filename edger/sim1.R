#!/usr/bin/env Rscript
# n the refman for edgeR these are only two mentions of simulation and this is one of them.
library(edgeR)
library(Cairo)

# Simulate DGE data for 1000 genes and 6 samples.
# Samples are in two groups
# First 200 genes are differentially expressed in second group
ngenes <- 1000
nlib <- 6 # aka. nsamps
counts <- matrix(rnbinom(ngenes*nlib, size=1/10, mu=20),ngenes,nlib)
rownames(counts) <- paste("gene",1:ngenes, sep=".")

# two groups which nehave differently
group <- gl(2,3,labels=c("Grp1","Grp2"))
counts[1:200,group=="Grp2"] <- counts[1:200,group=="Grp2"] + 10
y <- DGEList(counts,group=group)
y2 <- calcNormFactors(y)
# without labels, indexes of samples are plotted.
col <- as.numeric(group)

CairoPNG("sim1a.png", 800, 800)
mds <- plotMDS(y2, top=200, col=col)
show(mds)
dev.off()

# or labels can be provided, here group indicators:
CairoPNG("sim1b.png", 800, 800)
plotMDS(mds, col=col, labels=group)
dev.off()
