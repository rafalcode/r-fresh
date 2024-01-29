#!/usr/bin/env Rscript
# this script does what? Using mvtnorm package
library(mvtnorm)

# this is quite a good way to prefine a correlation between two random vairables
sigma <- matrix(c(4,2,2,3), ncol=2)
sigma <- matrix(c(7,2,2,3), ncol=2)
x <- rmvnorm(n=500, mean=c(1,2), sigma=sigma)

