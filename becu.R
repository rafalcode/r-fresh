#!/usr/bin/env Rscript
# generate beta distribution curves
library(ggplot2)
library(Cairo)

args <- commandArgs(trailingOnly = TRUE)
numargs <- length(args)
enumargs <- 2 # expected numargs
if(numargs != enumargs) {
    print("This script generates a beta distribution density curve, needing two arguments, shape1 and shape2")
    warning("Stopping right here")
}

sh1 <- as.numeric(args[1])
sh2 <- as.numeric(args[2])

sq <- seq(0, 1, by=0.05)
bvals <- dbeta(sq, sh1, sh2) # uniform beta vals (the x is uniform, not the y!)

CairoPNG("db0.png", 800, 800)
# lines connecting the dots (the type value is key here)
# plot(sq, bvals, pch=20, cex=1.5, col="firebrick", type="b")
plot(sq, bvals, pch=20, cex=1.5, col="firebrick")
lines(sq, bvals, col="grey")
# with lines you have more control over colour of the lines
dev.off()
