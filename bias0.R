#!/usr/bin/env Rscript
# this script does what? various ways of having bias
library(ggplot2)
library(Cairo)

# Set seed for reproducibility
set.seed(123)

# Number of data points
n <- 1000

# Parameters for the gamma distribution
shape <- 2   # Shape parameter (controls the shape of the distribution)
rate <- 0.5  # Rate parameter (controls the rate of decay)

# Generate data from a gamma distribution
biasda <- rgamma(n, shape = shape, rate = rate) # biased data
