#!/usr/bin/env Rscript
# instead of geom_bar() or geom_col() for the bars, 
# we actually use geom_segment() the approach is like points actually, but we just put 
# horizontal line segments.
library(tidyverse)
library(Cairo)

d <- data.frame(x=1:3, y = c(12, -9, 7))
rownames(d) <- sprintf("%02i", as.integer(rownames(d)))
lx <- dim(d)[1]
minl <- min(d$y)
barwidth <- 0.6
bd2 <- barwidth/2

# interpolate values from zero to y
vals <- lapply(d$y, function(y) seq(minl, y, by = 0.05))
y <- unlist(vals)
# create corresponding number of x values
mid <- rep(d$x, lengths(vals))
# so what we're going to do is generate horizontal line segmetns and so fill a bar
# stacking them upon each other to reach the required value.
d2 <- data.frame(x = mid - bd2, xend = mid + bd2, y = y, yend = y)

CairoPNG("fname.png", 800, 800)
ggplot(data = d2, aes(x = x, xend = xend, y = y, yend = yend, color = y)) +
  geom_segment(size = 1) +
  scale_color_gradient2(low = "red", mid = "white", high = "green", midpoint = 0)
dev.off()
