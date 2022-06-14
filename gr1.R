#!/usr/bin/env Rscript
# instead of geom_bar() or geom_col() for the bars, 
# we actually use geom_segment() the approach is like points actually, but we just put 
# horizontal line segments.
library(tidyverse)
library(Cairo)

d <- data.frame(x=1:3, y = c(12, -39, 7))
rownames(d) <- sprintf("%02i", as.integer(rownames(d)))
minl <- min(d$y)
maxl <- min(d$y)
if(abs(minl) > maxl) {
    ra <- abs(minl)
} else {
    ra <- maxl
}

barwidth <- 0.8
bd2 <- barwidth/2

# interpolate values from zero to y
vals <- lapply(d$y, function(y) seq(0, y, by = sign(y)*.025))
y <- unlist(vals)
# create corresponding number of x values
mid <- rep(d$x, lengths(vals))
# so what we're going to do is generate horizontal line segmetns and so fill a bar
# stacking them upon each other to reach the required value.
d2 <- data.frame(x = mid - bd2, xend = mid + bd2, y = y, yend = y)
# the y's and yend's are the same because the segmetn only moves on x axis,, its y coord is always the same
CairoPNG("fname.png", 800, 800)
ggplot(data = d2, aes(x = x, xend = xend, y = y, yend = yend, color = y)) +
  geom_segment(size = 1) +
  # scale_color_gradient2(low = "red", mid = "white", high = "green", midpoint = 0) +
  scale_color_gradient2(low = "midnightblue", mid = "grey", high = "green4", midpoint = 0) +
  theme(legend.position="none", panel.background = element_blank(), 
        axis.line.y=element_line(size = 0.5, linetype = "solid", colour = "black"),
        panel.grid.major = element_blank()) +
  ylim(-ra, ra) +
#   geom_hline(yintercept=0) +
  geom_rect(aes(xmin = 1, xmax = 3, ymin = 10, ymax = 15))
dev.off()
