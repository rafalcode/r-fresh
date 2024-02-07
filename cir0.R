#!/usr/bin/env Rscript
# this script does what? trying circos
library(circlize)
library(Cairo)

CairoPNG("cir0.png", 1400, 1400)
# circos.clear()
# circos.par$start.degree = 130 # yes must be before initialization
# circos.par$track.height=1.7
circos.initializeWithIdeogram(species="hg38", track.height=mm_h(4))
# circos.clear()
text(0, 0, "default", cex = 1)
dev.off()
