#!/usr/bin/env Rscript
# from https://jdblischak.github.io/dc-bioc-limma/batch-effects.html
# Blischak uses rembatchfx
# Bli
# this script does what?
# Note the problem with this script if you want to learn about removeBatchEffect, is that is conevrts to eset strucutre which
# makes it quite opaque, which ia  little disaapointing.
library(ggplot2)
library(Cairo)
library(dplyr)
library(limma)
library(edgeR)
# Have to load Biobase after dplyr so that exprs function works
library(Biobase)

# Accounting for technical batch effects: Response to infection with Mycobacterium tuberculosis
# John Blischak
# 2018-08-09

# An example of diagnosing and correcting batch effects from one of my own studies on the response to infection with Mycobacterium tuberculosis (paper, code, data).

# Download data.
# file_url <- "https://bitbucket.org/jdblischak/tb-data/raw/bc0f64292eb8c42a372b3a2d50e3d871c70c202e/counts_per_sample.txt"
# already downloaded:
fn <- "counts_per_sample.txt"
# full <- read.delim(file_url, stringsAsFactors = FALSE)
full <- read.delim(fn, stringsAsFactors = FALSE)

# dim(full)
# [1]   156 19419
# the first 5 or 6 columns are metadata, then on the header we have names of genes or features that's why there
# are so many columns
# the samples then, are in the columns

# Convert to ExpressionSet.
full <- full[order(full$dir), ]
rownames(full) <- paste(full$ind, full$bact, full$time, sep = ".")
x <- t(full[, grep("ENSG", colnames(full))])
p <- full %>% select(ind, bact, time, extr, rin)
stopifnot(colnames(x) == rownames(p))
eset <- ExpressionSet(assayData = x, phenoData = AnnotatedDataFrame(p))

# Filter lowly expressed genes.
keep <- rowSums(cpm(exprs(eset)) > 1) > 6
# sum(keep)
# [1] 12728

eset <- eset[keep, ]
# dim(eset)
# Features  Samples 
#    12728      156 

# Normalize with TMM.
norm_factors <- calcNormFactors(exprs(eset))
exprs(eset) <- cpm(exprs(eset), lib.size = colSums(exprs(eset)) * norm_factors, log = TRUE)
# plotDensities(eset, legend = FALSE)

# Clean up phenotype data frame to focus on early versus late timepoint for this example.
pData(eset)[, "infection"] <- ifelse(pData(eset)[, "bact"] == "none", "con", "inf")
pData(eset)[, "time"] <- ifelse(pData(eset)[, "time"] == 4, "early", "late")
pData(eset)[, "batch"] <- sprintf("b%02d", pData(eset)[, "extr"])

# table(pData(eset)[, c("time", "batch")])
#        batch
# time    b01 b02 b03 b04 b05 b06 b07 b08 b09 b10 b11 b12 b13
#   early   4   4   4   4   4   4   4   4   4   4   4   4   6
#   late    8   8   8   8   8   8   8   8   8   8   8   8   6

# Remove batch effect
# First Visualize principal components 1 and 2 for the original data.
CairoPNG("bef_rbe.png", 800, 800)
plotMDS(eset, labels = pData(eset)[, "time"], gene.selection = "common")
dev.off()
browser()
# Remove the effect of the technical variables: batch (discrete) and RIN (continuous; a measure of RNA quality).
exprs(eset) <- removeBatchEffect(eset, batch = pData(eset)[, "batch"], covariates = pData(eset)[, "rin"])

# Visualize principal components 1 and 2 for the corrected data.
CairoPNG("aft_rbe.png", 800, 800)
plotMDS(eset, labels = pData(eset)[, "time"], gene.selection = "common")
dev.off()
