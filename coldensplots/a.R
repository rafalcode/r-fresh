#!/usr/bin/env Rscript
# this script does what?
library(ggplot2)
library(Cairo)

CairoPNG("fname.png", 800, 800)
m <- ggplot(faithful, aes(x = eruptions, y = waiting)) +
 geom_point() +
  xlim(0.5, 6) +
   ylim(40, 110)
# contour lines
m + geom_density_2d()
dev.off()
