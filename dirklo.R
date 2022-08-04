#!/usr/bin/env Rscript
# Dirk's StOv post:
# https://stackoverflow.com/questions/3480388/how-to-fit-a-smooth-curve-to-my-data-in-r
# HOWEVER, this is not actually what I'm looking for,
# this (predict(loess()) is a standard model fitting exercise
# I want it to go through the points!
library(Cairo)

# I like loess() a lot for smoothing:

x <- 1:10
y <- c(2,4,6,8,7,12,14,16,18,20)
lo <- loess(y~x)

CairoPNG("dirk0.png", 800, 800)
plot(x,y)
lines(predict(lo), col='steelblue', lwd=2)
dev.off()

# Venables and Ripley's MASS book has an entire section on smoothing that also covers splines and polynomials -- but loess() is just about everybody's favourite.
