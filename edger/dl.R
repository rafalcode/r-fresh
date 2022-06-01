#!/usr/bin/env Rscript
# just the edgeR DGEList proof of concept, from the given example.
library(edgeR)
library(ggplot2)
library(Cairo)

ngenes <- 1000
nsamples <- 4

# we create a matrix of ngenes rows and nsamples column with random count values from rnbinom()
Counts <- matrix(rnbinom(ngenes*nsamples,mu=5,size=2),ngenes,nsamples)

rownames(Counts) <- 1:ngenes
y <- DGEList(counts=Counts, group=rep(1:2,each=2))
# dim(y)
# colnames(y)
browser()
y$samples
y$genes <- data.frame(Symbol=paste0("Gene",1:ngenes))

# Cairo image template
# CairoPNG("fname.png", 800, 800)
# # put plot command here
# dev.off()
