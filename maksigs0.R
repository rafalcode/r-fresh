#!/usr/bin/env Rscript
# this script does what? Maksi's geneset enrich
library(ggplot2)
library(Cairo)
library(here)
library(minfi)
library(paletteer)
library(limma)
library(reshape2)
library(missMethyl)
library(ggplot2)
library(glue)
library(tibble)
library(dplyr)
library(curatedTCGAData)
library(MultiAssayExperiment)
library(TCGAutils)
library(DMRcate)
library(patchwork)
source("maksiutility.R")

lv <- ls()
if(length(grep("^brca$", lv))!=1) {
    if(!file.exists("savedbrca.rds")) {
        brca <- curatedTCGAData(diseaseCode="BRCA", assays="Methylation_methyl450", version="2.1.1", dry.run=F)
    } else {
        brca <- readRDS("savedbrca.rds")
    }
}

patts <- list(normals = "TCGA-[A-Z0-9]{2}-[A-Z0-9]{4}-11")
hits <- makeHitList(brca, patts)
# brca <- MultiAssayExperiment::splitAssays(brca, c("11")) # extract only the normal samples
brca2 <- MultiAssayExperiment::splitAssays(brca, hits)
exp <- MultiAssayExperiment::experiments(brca2)[[1]]

meta <- colData(brca2)
betas <- as.matrix(assay(exp)) # should be 97 samples/columns .. these are the normals
colnames(betas) <- substr(colnames(betas), 1, 12)
m <- match(colnames(betas), meta$patientID)
meta <- meta[m, ]

betasNoNA <- betas[rowSums(is.na(betas)) == 0, ]
# unfiltered, unnecessary:
# mds <- plotMDS(betasNoNA, gene.selection = "common", plot = FALSE)
# dat <- tibble(x = mds$x, y = mds$y, gender = meta$gender)

# a bit of filtering
betasFlt <- rmSNPandCH(betasNoNA, rmXY = FALSE, rmcrosshyb = TRUE)
# dim(betasFlt) [1] 371789     97

# mds <- plotMDS(betasFlt, gene.selection = "common", plot = FALSE)
dat <- tibble(x = mds$x, y = mds$y, gender = meta$gender)

# logn forming:
dat <- as_tibble(melt(betasFlt))
colnames(dat) <- c("cpg", "ID", "beta")

numModes <- apply(betasFlt, 2, find_num_modes, adjust = 1.35)

dat <- as_tibble(melt(betasFlt))
colnames(dat) <- c("cpg", "ID", "beta")
dat$mode <- rep(as.character(numModes), each = nrow(betasFlt))

# head(meta[, 1:5])

# Cairo image template
# CairoPNG("fname.png", 800, 800)
# put plot command here
# dev.off()
