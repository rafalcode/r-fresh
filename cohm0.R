#!/usr/bin/env Rscript
# ComplexHeatmap tests, the anno_lines function.
# these are long thin recngular box to add onto heatmaps.
# i.e. as supplementary annotations
library(ComplexHeatmap)
library(Cairo)

anno = anno_lines(runif(10))
CairoPNG("cohm00.png", 800, 800)
draw(anno, test = "anno_lines")
dev.off()

anno = anno_lines(cbind(c(1:5, 1:5), c(5:1, 5:1)), gp = gpar(col = 2:3))
CairoPNG("cohm01.png", 800, 800)
draw(anno, test = "matrix")
dev.off()

anno = anno_lines(cbind(c(1:5, 1:5), c(5:1, 5:1)), gp = gpar(col = 2:3),
add_points = TRUE, pt_gp = gpar(col = 5:6), pch = c(1, 16))
CairoPNG("cohm02.png", 800, 800)
draw(anno, test = "matrix")
dev.off()
