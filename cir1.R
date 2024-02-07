#!/usr/bin/env Rscript
# this script does what? trying circos
library(circlize)
library(Cairo)

set.seed(123)

CairoPNG("cir1.png", 1400, 1400)
circos.initializeWithIdeogram(plotType = NULL)
circos.track(ylim = c(0, 1), panel.fun = function(x, y) {
    chr = CELL_META$sector.index
    xlim = CELL_META$xlim
    ylim = CELL_META$ylim
    circos.rect(xlim[1], 0, xlim[2], 1, col = rand_color(1))
    circos.text(mean(xlim), mean(ylim), chr, cex = 2.7, col = "white",
        facing = "inside", niceFacing = TRUE)
# }, track.height = 0.25, bg.border = NA)
}, track.height = mm_h(3), bg.border = NA)
text(0, 0, "default", cex = 1)
dev.off()
