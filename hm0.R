#!/usr/bin/env Rscript
# this script does what?
library(ggplot2)
library(Cairo)

# heatmaps can be as simple as:
x0 = data.matrix(UScitiesD, rownames.force = TRUE)
CairoPNG("fname0.png", 800, 800)
heatmap(x0, main = "Distances between US cities")
dev.off()
# as is, there a problem, perhaps its integers values?



ds = data.frame(rnorm(5, 50, 20),rnorm(5, 50, 20),rnorm(5, 50, 20),rnorm(5, 50, 20))
rn = c("Arm","Leg","Chest","Gut","Head")
cn = c("Ann","Bob","Tom","Joy")
x = data.matrix(ds, rownames.force = FALSE)

CairoPNG("fname2.png", 800, 800)
heatmap(x, labRow=rn, labCol=cn, main = "Test Heat Map")
dev.off()
