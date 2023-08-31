#!/usr/bin/env Rscript
# this script does what? Mixed models from chatgpt
library(mixtools)
library(Cairo)

# Number of data points
n <- 500

# True underlying means and standard deviations of the two normal distributions
true_mean1 <- 5
true_sd1 <- 1
true_mean2 <- 12
true_sd2 <- 2

# Generate synthetic data from two normal distributions
data1 <- rnorm(n, mean = true_mean1, sd = true_sd1)
data2 <- rnorm(n, mean = true_mean2, sd = true_sd2)

# Combine the data from both distributions
data <- c(data1, data2)

# Fit a two-component Gaussian mixture model to the data
fit <- normalmixEM(data, k = 2)

# Print the estimated parameters of the mixture model
# print(summary(fit))

# Plot the histogram of the data along with the fitted mixture model
# Cairo image template
CairoPNG("mx0.png", 800, 800)
plot(fit, which = 2)
dev.off()
