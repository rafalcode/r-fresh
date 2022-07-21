#!/usr/bin/env Rscript
# this script does what?
library(ggplot2)
library(Cairo)
library(Harman)
library(HarmanData)

data(OLF)
expt <- olf.info$Treatment
batch <- olf.info$Batch
olf.harman <- harman(olf.data, expt, batch)
CairoPNG("harp.png", 800, 800)
plot(olf.harman)
dev.off()
olf.data.corrected <- reconstructData(olf.harman)

## Reading from a csv file
datafile <- system.file("extdata", "NPM_data_first_1000_rows.csv.gz", package="Harman")
infofile <- system.file("extdata", "NPM_info.csv.gz", package="Harman")
datamatrix <- read.table(datafile, header=TRUE, sep=",", row.names="probeID")
batches <- read.table(infofile, header=TRUE, sep=",", row.names="Sample")
res <- harman(datamatrix, expt=batches$Treatment, batch=batches$Batch)

CairoPNG("hara.png", 800, 800)
arrowPlot(res, 1, 3)
dev.off()
