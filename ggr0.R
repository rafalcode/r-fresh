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

Cairo(1000, 500, "ggr0.png")

dat <- subset(mtcars, wt > 2.75 & wt < 3.45)
dat$car <- rownames(dat)

p <- ggplot(dat, aes(wt, mpg, label = car)) + geom_point(color = "red")
p1 <- p + geom_text() + labs(title = "geom_text()")
p2 <- p + geom_text_repel() + labs(title = "geom_text_repel()")

gridExtra::grid.arrange(p1, p2, ncol = 2)
dev.off()
