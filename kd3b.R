#!/usr/bin/env Rscript
# this script does what?
# looking to replace https://commons.wikimedia.org/wiki/File:Comparison_of_1D_histogram_and_KDE.png
library(Cairo)
set.seed(2)

x <- seq(-8, 8, .2)
rx <- c(-2.1, -1.3, -.4, 1.9, 5.1, 6.2) # our received values - only six! when smoothened distri is likely to be be odd.

d <- density(rx, bw=1.3) # eyeballed bandwidth

