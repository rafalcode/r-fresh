#!/usr/bin/env Rscript
# this script does what?
library(ggplot2)
library(Cairo)
library(minfi)
library(mCSEA)
library(FlowSorted.Blood.450k)

# Cairo image template
# CairoPNG("fname.png", 800, 800)
# put plot command here
# dev.off()

minfiDataDir <- system.file("extdata", package = "minfiData")
targets <- read.metharray.sheet(minfiDataDir, verbose = FALSE)
RGset <- read.metharray.exp(targets = targets)

data(mcseadata)
cellCounts = estimateCellCounts(RGset) # a minfi function.
myRank <- rankProbes(betaTest, phenoTest, refGroup = "Control")
