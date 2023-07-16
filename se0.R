#!/usr/bin/env Rscript
# this script does what? Exercises in SummarizedExperiemnt package, much beloved by nfcore crowd.
library(SummarizedExperiment)
library(Cairo)

# I think we can say that SummarizedExperiemtn (aka. "se") packages thec ounts data, the phenotypical (covatiate) data
# and what probably is the annotation data altogether,

nrows <- 20; ncols <- 6

counts <- matrix(runif(nrows * ncols, 1, 1e4), nrows)
# floating point counts to be sure

rowRanges <- GRanges(rep(c("chr1", "chr2"), c(5, 15)),
                     IRanges(sample(1000L, 20), width=100),
                     strand=Rle(c("+", "-"), c(12, 8)),
                     seqlengths=c(chr1=1800, chr2=1300))

colData <- DataFrame(Treatment=rep(c("ChIP", "Input"), 3),
                     row.names=LETTERS[1:6])

# so now we have the three things to build a summarized experiement.
rse <- SummarizedExperiment(assays=SimpleList(counts=counts),
                            rowRanges=rowRanges, colData=colData)
cvg <- coverage(rse)

# Cairo image template
# CairoPNG("fname.png", 800, 800)
# put plot command here
#dev.off()
