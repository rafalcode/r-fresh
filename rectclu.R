#!/usr/bin/env Rscript
# this script does what?
library(Cairo)

hca <- hclust(dist(USArrests))

CairoPNG("fname.png", 800, 800)
plot(hca)
rect.hclust(hca, k = 3, border = "red")
# x <- rect.hclust(hca, h = 50, which = c(2,7), border = 3:4)
dev.off()
