#!/usr/bin/env Rscript
library(affy)
library(gplots)
library(Cairo)

data(SpikeIn)
pms <- SpikeIn@pm

# just the data, scaled across rows
CairoPNG("aff.png", 800, 800)
heatmap.2(pms, col=rev(heat.colors(16)), main="SpikeIn@pm", xlab="Relative Concentration", ylab="Probeset", scale="row")
dev.off()

# fold change vs "12.50" sample
data <- pms / pms[, "12.50"]
data <- ifelse(data>1, data, -1/data)
CairoPNG("aff2.png", 800, 800)
heatmap.2(data, breaks=16, col=redgreen, tracecol="blue", main="SpikeIn@pm Fold Changes\nrelative to 12.50 sample", xlab="Relative Concentration", ylab="Probeset")
dev.off()
