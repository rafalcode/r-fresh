#!/usr/bin/env Rscript
# this script does what? Construct a SummarizedExpermient object.
library(SummarizedExperiment)

nrows <- 200
ncols <- 6
counts <- matrix(runif(nrows * ncols, 1, 1e4), nrows)
cnamsq <- rep(c("chr1", "chr2"), c(50, 150)) # chromosome name sequence, length 200
runifs <- floor(runif(200, 1e5, 1e6)) # l. 200
rastra <- sample(c("+", "-"), 200, TRUE) # l. 200 random strands
rowRanges <- GRanges(cnamsq, IRanges(runifs, width=100), strand=rastra, feature_id=sprintf("ID%03d", 1:200))

# colData <- DataFrame(Treatment=rep(c("ChIP", "Input"), 3), row.names=LETTERS[1:6])
colData <- DataFrame(Tissue=rep(c("Norm", "Tumour"), 3), row.names=LETTERS[1:6])

se1 <- SummarizedExperiment(assays=list(counts=counts), rowRanges=rowRanges, colData=colData)
## class: RangedSummarizedExperiment 
## dim: 200 6 
## metadata(0):
## assays(1): counts
## rownames: NULL
## rowData names(1): feature_id
## colnames(6): A B ... E F
## colData names(1): Treatment

# A SummarizedExperiment can be constructed with or without supplying a DataFrame for the rowData argument:
se2 <- SummarizedExperiment(assays=list(counts=counts), colData=colData)
## class: SummarizedExperiment 
## dim: 200 6 
## metadata(0):
## assays(1): counts
## rownames: NULL
## rowData names(0):
## colnames(6): A B ... E F
## colData names(1): Treatment
