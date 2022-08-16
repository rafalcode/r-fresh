#!/usr/bin/env Rscript
# I want to tyr Summarized Exeriments objects.
library(SummarizedExperiment)
nrows <- 200
ncols <- 6
randcounts <- matrix(runif(nrows * ncols, 1, 1e4), nrows)

chrnams <- rep(c("chr1", "chr2"), c(nrows/4, 3*nrows/4))
nums <- floor(runif(200, 1e5, 1e6))
irnums <- IRanges(nums, width=100)
rowRanges <- GRanges(chrnams, irnums,
                     strand=sample(c("+", "-"), 200, TRUE),
                     feature_id=sprintf("ID%03d", 1:nrows))
colData <- DataFrame(Treatment=rep(c("ChIP", "Input"), 3),
                     row.names=LETTERS[1:6])

summex <- SummarizedExperiment(assays=list(counts=randcounts), rowRanges=rowRanges, colData=colData)
