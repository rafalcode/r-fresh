#!/usr/bin/env Rscript
# this script does what?
library(ggplot2)
library(Cairo)
library(ggnewscale)

# Cairo image template

# Equivalent to melt(volcano), but we don't want to depend on reshape2
topography <- expand.grid(x = 1:nrow(volcano), y = 1:ncol(volcano))
topography$z <- c(volcano)

# point measurements of something at a few locations
measurements <- data.frame(x = runif(30, 1, 80), y = runif(30, 1, 60), thing = rnorm(30))

CairoPNG("fname.png", 800, 800)
ggplot(mapping = aes(x, y)) +
      geom_contour(data = topography, aes(z = z, color = stat(level))) +
      # Color scale for topography
      scale_color_viridis_c(option = "D") +
      # geoms below will use another color scale
      new_scale_color() +
      geom_point(data = measurements, size = 3, aes(color = thing)) +
      # Color scale applied to geoms added after new_scale_color()
      scale_color_viridis_c(option = "A")
# put plot command here
dev.off()
