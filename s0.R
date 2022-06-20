#!/usr/bin/env Rscript
# this script does what?
library(ggplot2)
library(Cairo)
# Cairo image template
CairoPNG("fname.png", 800, 800)
p1 <- ggplot(mpg, aes(displ, hwy)) + geom_point() +
    # scale_y_continuous(label = c("two", "four", "six"))
    # scale_y_continuous(breaks = c(15, 20, 25, 30, 35, 40), label = c("a", "b", "c", "a", "b", "c"))
    scale_y_continuous(breaks = c(15, 20, 25, 30, 35, 40), label = c(2,1,0,0,1,2))
show(p1)
dev.off()
