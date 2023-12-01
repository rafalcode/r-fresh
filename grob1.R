#!/usr/bin/env Rscript
# I'm trying to understand this grob thing and how to apply it
library(Cairo)
library(grid)
library(ggplot2)

CairoPNG("grob1.png", 800, 800)

d <- ggplot(diamonds, aes(carat, price, fill = ..density..)) +
      xlim(0, 2) +
      stat_binhex(na.rm = TRUE) +
      # opts(aspect.ratio = 1) +
      facet_wrap(~ color, ncol = 4)
show(d)
dev.off()
