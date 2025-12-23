#!/usr/bin/env Rscript
# this script does what?
library(Cairo)
set.seed(2)

svgpts <- function() {
    # these are the points gleaned from the SVG
    pts <- c(71.922, 71.922, 78.321, 90.281, 93.718, 94.994, 95.73, 98.748, 100.927, 101.365, 102.026, 105.678, 108.229, 109.096, 109.124, 109.618, 110.177, 111.416, 111.77, 113.726, 116.836, 117.498, 118.522, 119.137, 120.292, 120.86, 121.81, 122.08, 122.955, 123.431, 124.707, 125.656, 127.669, 128.814, 130.463, 132.308, 132.484, 132.689, 133.183, 133.844, 134.123, 135.865, 136.126, 136.703, 137.076, 137.645, 137.914, 138.743, 141.752, 143.065, 143.494, 143.988, 144.351, 144.872, 145.096, 145.813, 147.2, 148.496, 148.831, 149.11, 150.675, 151.495, 152.91, 155.091, 155.397, 155.528, 156.906, 156.935, 157.018, 157.986, 159.049, 164.014, 164.209, 165.066, 165.812, 167.013, 167.683, 169.089, 169.247, 169.956, 170.216, 171.754, 171.968, 172.331, 172.387, 175.545, 175.666, 176.765, 176.895, 177.239, 178.386, 181.552, 183.126, 183.825, 187.709, 192.98, 194.304, 197.023, 204.475, 212.02)
    
    # first is repeated, leavin gus only 99 ... must make do.
    pts <- pts[2:100]
    
    # Eye up zero element (48) and set first value to be around 2.3
    dv <- pts[1]/2.4
    p2 <- (pts-pts[48])/dv
}

x <- seq(-3.1,3.1,.1)
# rx <- rnorm(100)
rx <- svgpts()
stroke <-c("#333399", "#009900", "#FF6633", "#CCCCCC")

# CairoPNG("kd2.png", 800, 800)
CairoSVG("kd2.svg", width=10, height=8)
# par(xaxs="i")

plot(x, dnorm(x), "l", ylim=c(0,.5), col=stroke[4], lwd=5, ylab="Density") # that "l" is for lines
lines(density(rx,bw=0.3), col=stroke[1], lwd=3.0)
lines(density(rx,bw=0.1), col=stroke[2], lwd=3.0)
lines(density(rx,bw=0.05), col=stroke[3], lwd=3.0)
legend(-3.1, .5, legend=c("Reference", "0.3", "0.1", "0.05"), col=stroke[4:1], lty=1, cex=1.6, bty="n", lwd=c(5, 3, 3, 3))
# bty is border/box type I reckon. N means no box.
rug(rx)
# Cairo image template
# put plot command here
dev.off()
