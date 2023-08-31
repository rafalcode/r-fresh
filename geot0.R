#!/usr/bin/env Rscript
# this script does what?
library(GeoTcgaData) # for the differential_RNA() func.
library(SummarizedExperiment)
library(Cairo)

# use user-defined data
nrows <- 25
ncols <- 16
df <- matrix(rnbinom(nrows*ncols, mu = 4, size = 10), nrows, ncols)
df <- as.data.frame(df)
rownames(df) <- paste0("G", 1:nrows)
colnames(df) <- paste0("S", 1:ncols)
group <- sample(c("group1", "group2"), ncols, replace = T) # unfort, randomly assigned.
res0 <- differential_RNA(counts = df, group = group, filter = FALSE, method = "Wilcoxon")

# use SummarizedExperiment object input
df2 <- matrix(rnbinom(nrows*ncols, mu = 4, size = 10), nrows, ncols) # again, like.
rownames(df2) <- paste0("G", 1:nrows)
colnames(df2) <- paste0("S", 1:ncols)
group2 <- sample(c("group1", "group2"), 16, replace = TRUE)

colData <- S4Vectors::DataFrame(row.names = paste0("S", 1:ncols), group = group2)
sedata <- SummarizedExperiment::SummarizedExperiment(
         assays=S4Vectors::SimpleList(counts=df2),
         colData = colData)
res1 <- differential_RNA(counts = sedata, groupCol = "group", filter=F, method="Wilcoxon") 
res2 <- differential_RNA(counts = sedata, groupCol = "group", filter=F, method="DESeq2") 

# meth oriented example. no action on it though.
nrows <- 200
ncols <- 20
counts <- matrix(
    runif(nrows * ncols, 1, 1e4), nrows,
    dimnames = list(paste0("cg",1:200),paste0("S",1:20)))


# Cairo image template
# CairoPNG("fname.png", 800, 800)
# put plot command here
# dev.off()
