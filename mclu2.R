#!/usr/bin/env Rscript
# used chag to get a simualtion code for the old faitful data set
library(ggplot2)
library(Cairo)
library(mclust)

gmm_model <- Mclust(faithful$eruptions, G=2)  # Fit mixture of 2 Gaussians
summary(gmm_model)  # Get estimated means, variances, and probabilities

# Simulate using the fitted model parameters
simerups_gmm <- rnorm(n, mean=gmm_model$parameters$mean, sd=sqrt(gmm_model$parameters$variance$sigmasq))

# Cairo image template
CairoPNG("mclu20.png", 800, 800)
# Plot comparison
par(mfrow=c(1,2))
hist(faithful$eruptions, breaks=20, probability=TRUE, main="Real Data", col="lightblue")
hist(simerups_gmm, breaks=20, probability=TRUE, main="Simulated Data (GMM)", col="lightgreen")


# put plot command here
dev.off()
