#!/usr/bin/env Rscript
# this script does what?
library(ggplot2)
library(Cairo)

diam <- diamonds
# dim is 53940    10
sa <- sample(1:nrow(diam), 200)

# Cairo image template
CairoPNG("diam.png", 800, 800)
ggp <- ggplot(diam[sa,], aes(carat)) +
     geom_density()
show(ggp)
dev.off()
