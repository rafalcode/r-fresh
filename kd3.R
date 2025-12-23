#!/usr/bin/env Rscript
# this script does what?
library(Cairo)
set.seed(2)

x <- seq(-8, 8, .2)
rx <- c(-2.1, -1.3, -.4, 1.9, 5.1, 6.2)

h <- hist(rx, breaks=4, plot=F)
d <- density(rx)

CairoPNG("kd3.png", 800, 800)
plot(h, border="blue", freq=F)
box()
rug(rx, lwd=3, ticksize=0.02)
# lines(density(rx,bw=0.3))
lines(density(d))
# lines(density(rx,bw=

# Cairo image template
# put plot command here
dev.off()
