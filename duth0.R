#!/usr/bin/env Rscript
# this script does what?
# exploring Brad Duthie stuff
# to wit: https://stirlingcodingclub.github.io/simulating_data/index.html
library(ggplot2)
library(Cairo)


# this is your typical correlation sinulation. Firs we start with the predictor variable x1
N   <- 10000
rho <- 0.8
# x1  <- rnorm(n = N, mean = 0, sd = 1)
x1  <- rnorm(N) # this is the same, as above are defaults.

# OK the depedent variable, this is a stdanrd way to do it.
# x2  <- (rho * x1) + sqrt(1 - rho*rho) * rnorm(n = N, mean = 0, sd = 1);
# amounts to:
rho2 <- rho*rho # come sto .09 in this case
sq1mr <- sqrt(1-rho2) # square root of 1 minus rho squared, comes to .95 in this case
# x2  <- (rho * x1) + sqrt(1 - rho*rho) * rnorm(N)
x2  <- (rho * x1) + sq1mr*rnorm(N)

# what you get here is an x1 of .3 suprise surprise with incredibly low p-value.
# of course you could say, of course, the normal distri "plays into its hands"
# HOWEVER, note how small the rsquared is!
# you can raise it, by bringing rho up to .8, this will give you an rqs of .635 (not that adjusted is the exact same!) 

# Cairo image template
# CairoPNG("fname.png", 800, 800)
# put plot command here
# dev.off()
