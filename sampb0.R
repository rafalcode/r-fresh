#!/usr/bin/env Rscript
# this script does what? this uses chatgpt to introduce sampling bias
library(ggplot2)
library(Cairo)

# Set the seed for reproducibility
# set.seed(123)

origsz <- 4000
subsz <- 40

# Generate a large list of numbers
origset <- rnorm(origsz, mean = 50, sd = 10)

# Define the bias factors
biased_factor <- 4  # Higher values mean higher probability of selection

# Calculate probs for biased sampling
probs <- (origset - min(origset))^biased_factor
probs <- probs / sum(probs)

# Sample the biased subset
biasedsub <- sample(origset, size = subsz, replace = FALSE, prob = probs)

# Print summary statistics
cat("Original Mean:", mean(origset), "\n")
cat("Biased Subset Mean:", mean(biasedsub), "\n")
cat("Original sd:", sd(origset), "\n")
cat("Biased Subset sd:", sd(biasedsub), "\n")

# Cairo image template
# CairoPNG("fname.png", 800, 800)
# put plot command here
# dev.off()
