#!/usr/bin/env Rscript
# this script s an exercise in basing one random variable 
# off another using a strength paramter.
# the strength parameter is the influence of the first random variable on the second
# the "rest" of the strength (1-sqrt(1-strength^2)) is intensity of the random aspect 
# of the second variable which is not dependent on the first.
# library(ggplot2)
library(Cairo)
N <- 50
nbreaks <- N/5

rr <- rnorm(N)
strength <- 0.95
relement <- 1-sqrt(1-strength^2)
rfacs <- relement*rnorm(N)
r2 <- strength*rr + rfacs

# A simple plot() will just give you index point on x and corresponding vector value on y, not so useful
# CairoPNG("rap0.png", 800, 800)
# plot(rr)
# dev.off()


CairoPNG("rap0.png", 800, 800)
hist(rr, breaks =nbreaks)
dev.off()

CairoPNG("rap1.png", 800, 800)
hist(rfacs, breaks =nbreaks)
dev.off()

CairoPNG("rap2.png", 800, 800)
hist(r2, breaks =nbreaks)
dev.off()

CairoPNG("rap3.png", 800, 800)
plot(x=rr, y=r2)
dev.off()
cat(paste0("Cor(rr,rfacs) is ", cor(rr,rfacs)), "\n")
cat(paste0("Cor(rr,r2) is ", cor(rr,r2)), "\n")


# Notes what I'm finding is that when strength is low, there is high correlation between rr and r2.
# I have to "raise" the strength to 0.95 to get a still correlation of 0.85

# So this method is not too great, and works in the opposite way I would have thought.
