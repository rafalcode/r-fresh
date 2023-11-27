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
cat(paste0("We take an rnorm randvar vector of length ", N, " and derive a new vector from it by adding a weakened (", strength, ") version of original to a new randvec multiplied by ", relement, ".\n"))
# I submitted this to ChatGPT and what was returned was reasonable:

# It looks like you've provided some R code that generates random numbers using the rnorm function and combines them to create a new random variable r2 with a specified correlation (strength) between rr and rfacs.
# Here's a brief breakdown of what this code does:
# rr <- rnorm(N): Generates N random numbers from a standard normal distribution and assigns them to the variable rr.
# strength <- 0.95: Specifies a correlation strength of 0.95 between rr and rfacs.
# relement <- 1-sqrt(1-strength^2): Calculates the scaling factor relement to achieve the desired correlation. This formula is based on the mathematical relationship between the correlation coefficient and the scaling factor for two variables.
# rfacs <- relement*rnorm(N): Generates N random numbers from a standard normal distribution, multiplies them by the relement scaling factor, and assigns the result to the variable rfacs. These values are intended to be used to introduce the desired correlation with rr.
# r2 <- strength*rr + rfacs: Combines rr and rfacs to create a new random variable r2 with the specified correlation strength.
# The code appears to be correctly implementing the concept of generating correlated random variables using the correlation coefficient strength. If you have any specific questions or need further clarification about this code, please feel free to ask.


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

CairoPDF("rap2.pdf", 800, 800)
hist(r2, breaks =nbreaks)
dev.off()

CairoPNG("rap3.png", 800, 800)
plot(x=rr, y=r2)
dev.off()
cat(paste0("Cor(rr,rfacs) is ", cor(rr,rfacs)), "\n")
cat(paste0("Cor(rr,r2) is ", cor(rr,r2)), "\n")

cat(paste0("Orig vec = ", paste0(rr, collapse=" "), "\n"))
cat(paste0("Derived vec = ", paste0(r2, collapse=" "), "\n"))

# Notes what I'm finding is that when strength is low, there is high correlation between rr and r2.
# I have to "raise" the strength to 0.95 to get a still correlation of 0.85

# So this method is not too great, and works in the opposite way I would have thought.
