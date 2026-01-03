#!/usr/bin/env Rscript
# this script does what?
# looking to replace https://commons.wikimedia.org/wiki/File:Comparison_of_1D_histogram_and_KDE.png
library(Cairo)
set.seed(2)

x <- seq(-8, 8, .2)
rx <- c(-2.1, -1.3, -.4, 1.9, 5.1, 6.2) # eyeballed values 

h <- hist(rx, breaks=4, plot=F)
d <- density(rx, bw=1.3) # eyeballed bandwidth

# CairoPNG("kd3.png", 1600, 800)
CairoSVG("kd3.svg", 16, 8)
par(mfrow=c(1,2))
par(lwd=2) # will both axes and the bar borders to this width by default.
# par(mar = c(5, 8, 4, 2))
par(mar = c(5, 5, 4, 2))
# plot(h, border="blue", lwd=2, col="white", freq=F, xlab="x", ylab="Density function", main="", cex.lab=1.5, xlim=c(-5.5, 11), ylim=c(0,.17)) # watch that freq option.
plot(h, border="blue", lwd=2, col="white", freq=F, xlab="x", ylab="Density function", main="", cex.lab=1.5, xlim=c(-5.5, 11), ylim=c(0,.17)) # watch that freq option.

box() # make a box in stead of isolated x and y
rug(rx, lwd=3, ticksize=0.03)
# lines(density(rx,bw=0.3))
# lines(d)
# plot(d, type="l", lwd=2, col="blue", xlab="x", ylab="Density function", main="", cex.lab=1.5, xlim=c(-5.5, 11), ylim=c(0,.17))
plot(d, type="l", col="blue", xlab="x", ylab="Density function", main="", cex.lab=1.5, xlim=c(-5.5, 11), ylim=c(0,.17))
# now draw a normal above each data point, to illustrate KDE.

# all possible points.
ss0 <- seq(rx[1]-.5,rx[6]+.5,.1)
m <- matrix(rep(0, length(ss0)*length(rx)), nrow=length(ss0), ncol=length(rx)) ## on which you will perform rowSums

for(i in 1:length(rx)) {
    w <- which(x>(rx[i]-2.5) & x<(rx[i]+2.5)) # find the x values of interest to this particular data point.
    x0 <- x[w]
    d0 <- .1*dnorm(x0, mean=rx[i], sd=1)
    m[w,i] <- d0 # collect in matrix
    lines(x0, y=d0, type="l", lwd=2, lty=2, col="red") # 5 is long dash
}
lines(ss0, y=rowSums(m), type="l", lwd=2, lty=2, col="green") # 5 is long dash
rug(rx, lwd=3, ticksize=0.03)
# lines(density(rx,bw=

# You should be able to plot rowSums(m) now and it should look like d.
# Cairo image template
# put plot command here
dev.off()
