#!/usr/bin/env Rscript
# this script does what? trying circos
library(circlize)
library(Cairo)

set.seed(123)

human_cytoband = read.cytoband(species = "hg38")$df
human_cytoband[ ,1] = paste0("human_", human_cytoband[, 1])

chromosome.index = paste0("human_chr", c(1:22, "X", "Y"))

CairoPNG("cir2.png", 1400, 1400)
circos.initializeWithIdeogram(human_cytoband, plotType = NULL, chromosome.index = chromosome.index)
circos.track(ylim = c(0, 1), panel.fun = function(x, y) {
    circos.text(CELL_META$xcenter, CELL_META$ylim[2] + mm_h(2)
        gsub(".*chr", "", CELL_META$sector.index), cex=.6, niceFacing = TRUE)
}, track.height = .5, cell.padding = c(0, 0, 0, 0), bg.border = NA)

# circos.genomicIdeogram(human_cytoband, track.height=.1)
circos.genomicIdeogram(human_cytoband)
dev.off()
