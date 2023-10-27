#!/usr/bin/env Rscript
# this script does what?
# this is the example in the WGCNA refman which explains how WGCNA does it faster.
library(WGCNA)
library(DGCA)

set.seed(10)
nrow <- 100
ncol <- 500

data <- matrix(rnorm(nrow*ncol), nrow, ncol);
colnames(data) <- paste0("S", 1:ncol)
rownames(data) <- paste0("G", 1:nrow)

## First test: no missing data
st0 <- system.time( {corStd <- stats::cor(data)} )
# st1 <- system.time( {corFast <- WGCNA::cor(data)} )
corFast2 <- WGCNA::cor(data, use = "pairwise.complete.obs")
mcsig <- DGCA::matCorSig(corFast2, ncol)
all.equal(corStd, corFast)

# Here R's standard correlation performs very well.
# We now add a few missing entries.
data[sample(nrow, 10), 1] <- NA

# And test the correlations again...
st2 <- system.time( {corStd <- stats::cor(data, use='p')} );
st3 <- system.time( {corFast <- WGCNA::cor(data, use='p')} );
all.equal(corStd, corFast)
# Here the R's standard correlation slows down considerably
# while corFast still retains it speed. Choosing
# higher ncol above will make the difference more pronounced.
