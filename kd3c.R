#!/usr/bin/env Rscript
# this script does what?
# I want to try splicefun, see how different of close it is to density plot IN THE NAIVE MANNER.
# In the naive manner? well, you see density is based on kernels, little normal distributions places over each point.
# so splinefun needs more if it's going to get closer to density ... and may even need to use the desnity function itself first.
library(Cairo)
set.seed(2)

rx <- c(-2.1, -1.3, -.4, 1.9, 5.1, 6.2) # eyeballed values 

h <- hist(rx, breaks=4, plot=F)
d <- density(rx, bw=1.3) # eyeballed bandwidth
a <- approx(h$mids, h$density)
spl <- spline(d$x, d$y)

CairoPNG("kd3c.png", 1600, 800)
par(mfrow=c(1,2))
par(lwd=2) # will both axes and the bar borders to this width by default.
par(mar = c(5, 5, 4, 2))
# plot(d, type="l", col="blue", xlab="x", ylab="Density function", main="", cex.lab=1.5, xlim=c(-5.5, 11), ylim=c(0,.17))
plot(h, col="orange", xlab="x", ylab="hist() function", main="", cex.lab=1.5, xlim=c(-5.5, 11), ylim=c(0,.17))
lines(a, col="purple")
plot(a, type="l", lwd=2, col="purple", xlab="x", ylab="Approx of mids and h$dens", main="", cex.lab=1.5, xlim=c(-5.5, 11), ylim=c(0,.17)) # watch that freq option.
box() # make a box in stead of isolated x and y
rug(rx, lwd=3, ticksize=0.03)

plot(d, type="l", col="blue", xlab="x", ylab="Density function", main="", cex.lab=1.5, xlim=c(-5.5, 11), ylim=c(0,.17))
box() # make a box in stead of isolated x and y
rug(rx, lwd=3, ticksize=0.03)

dev.off()
