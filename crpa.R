#!/usr/bin/env Rscript
# this script does what?
library(ggplot2)
library(Cairo)

paletteFunc <- colorRampPalette(c("#0101ed", "#fbfb02"))
palette  <- paletteFunc(32);

# Cairo image template
CairoPNG("fname.png", 800, 800)
barplot(1:8, col=palette)
dev.off()
