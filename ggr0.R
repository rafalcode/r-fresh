#!/usr/bin/env Rscript
library(Cairo)
library(ggrepel)
# example from ggrepel's examples: Align labels on the top or bottom edge

# Use hjust to justify the text neatly:
# 
#     hjust = 0 for left-align
#     hjust = 0.5 for center
#     hjust = 1 for right-align
# 
# Sometimes the labels do not align perfectly. Try using direction = "x" to limit label movement to the x-axis (left and right) or direction = "y" to limit movement to the y-axis (up and down). The default is direction = "both".
# 
# Also try using xlim() and ylim() to increase the size of the plotting area so all of the labels fit comfortably.

set.seed(42)

Cairo(800, 800, "ggr0.png")
# ggplot(mtcars, aes(x = wt, y = 1, label = rownames(mtcars))) +
ggplot(mtcars, aes(x = 1, y = wt, label = rownames(mtcars))) +
  geom_point(color = "red") +
  geom_text_repel(
    force_pull   = 0, # do not pull toward data points
    nudge_y      = 0.05,
    direction    = "x",
    angle        = 90,
    hjust        = 0,
    segment.size = 0.2,
    max.iter = 1e4, max.time = 1) +
  xlim(1, 6) + ylim(1, 0.8) +
  theme(
    axis.line.y  = element_blank(),
    axis.ticks.y = element_blank(),
    axis.text.y  = element_blank(),
    axis.title.y = element_blank())
dev.off()
