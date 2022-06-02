#!/usr/bin/env Rscript
# just the edgeR DGEList proof of concept, from the given example.
library(edgeR)
library(ggplot2)
library(Cairo)

ndegenes <- 200
nngenes <- 800
ngenes <- nngenes+ndegenes
nsamples <- 6

# we create a matrix of ngenes rows and nsamples column with random count values from rnbinom()
Counts <- matrix(rnbinom(ngenes*nsamples,mu=5,size=2), ngenes,nsamples)
group <- rep(1:2,each=nsamples/2)
# Counts[nngenes+1:ngenes,2] <- rnbinom(ndegenes,mu=20,size=2)
# Counts[nngenes+1:ngenes,4] <- rnbinom(ndegenes,mu=20,size=2)
# Counts[nngenes+1:ngenes,6] <- rnbinom(ndegenes,mu=20,size=2)
Counts[(nngenes+1):ngenes,c(2,4,6)] <- rnbinom(3*ndegenes,mu=20,size=2)

rownames(Counts) <- 1:ngenes
y <- DGEList(counts=Counts, group=rep(1:2,each=nsamples/2))

#y, or the DGEList object has a norm/factors member but they are all at 1
# note they are "per column" so 6 of them here.
ynf <- calcNormFactors(y)
# now ynf's norm.factors are not 1 but show real weighting for each sample.
# seems clear this is for sequencing depth.
CPMs <- cpm(ynf, log=FALSE, prior.count=1)
logCPMs <- cpm(ynf, log=TRUE, prior.count=1)
