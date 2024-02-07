#!/usr/bin/env Rscript
# this script does what? trying circos this is from:
# https://jokergoo.github.io/circlize_book/book/initialize-genomic-plot.html#initialize-cytoband
library(circlize)
library(Cairo)

set.seed(123)

human_cytoband = read.cytoband(species = "hg38")$df
human_cytoband[ ,1] = paste0("human_", human_cytoband[, 1])

chromosome.index = paste0("human_chr", c(1:22, "X", "Y"))

CairoPNG("cir2.png", 600, 600)
circos.initializeWithIdeogram(human_cytoband, plotType = NULL, chromosome.index = chromosome.index)
circos.track(ylim = c(0, 1), panel.fun = function(x, y) {
    # circos.text(CELL_META$xcenter, CELL_META$ylim[2] + mm_y(2),
    circos.text(CELL_META$xcenter, CELL_META$ylim[2],
        gsub(".*chr", "", CELL_META$sector.index), cex=1.2, niceFacing = TRUE)
# }, track.height=mm_h(3), cell.padding = c(0, 0, 0, 0), bg.border = NA)
}, track.height=.025, cell.padding = c(0, 0, 0, 0), bg.border = NA)
# what are mm_y() and mm_h() ... niot mus musculus, but millimeters, h for horiz, y for y-axis?
# including mm_y(2) causes problem, possibly due to modification

circos.genomicIdeogram(human_cytoband, track.height=.075)
# circos.genomicIdeogram(human_cytoband)
dev.off()
