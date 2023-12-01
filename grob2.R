#!/usr/bin/env Rscript
# this script does what?
library(ggplot2)
library(Cairo)
library(grid)

# Cairo image template
CairoPNG("grob2.png", 800, 800)

my_circle <- circleGrob(name = "my_circle", x = 0.5, y = 0.5, r = 0.5, gp = gpar(col = "gray80", lty = 3))
# lty: line type I bet
grid.draw(my_circle)

my_rect <- rectGrob(x = 0.5, y = 0.5, width = 0.8, height = 0.3)
grid.draw(my_rect)
dev.off()
