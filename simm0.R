#!/usr/bin/env Rscript
# this s from Peters 2015 supplementary paper, but he doesn't give the ocde just a garbled version.
library(ggplot2)
library(Cairo)

# beta2 = beta1 + 0.20
# This means that
# beta2 ∼ Uniform(0.21, 0.99)
# For “up” regions, the base methylation level for the control samples was
# set as beta1, and for treatment samples as beta2. For “down” regions this
# allocation was reversed.
# 3. Structure of simulated data sets
# We simulated 10 replicate samples in each of the two sample groups, con-
#     trol and treatment. Thus each simulated data set was generated as a
# 2D array with 485,512 rows and 20 columns. The first 10 columns were
# control samples, the last 10 columns wer treatment samples.
# 4. Simulation of beta values: Inside DMRs
# First we simulate beta values in each of the DMRs selected in stage (2)
# above, using the two beta level values, one for control samples and one
# for treatment samples. Random values are generated for all samples and
# all probes in the region. A beta distribution is used to generate random
# data (function rbeta() in R). The beta distribution has two parameters,
# a and b. The mode of the beta(a, b) distribution is

# For inside a DMR
a <- 18 # mine!
b <- 80 # mine!
mu <- (a - 1) / (a + b - 2) # also called mode

# We're choosing parameters a and b so that this mode is equal to the specified beta ’level’ and so that
# a + b + 2 = K = 100
K <- 100
# This value was chosen to give a realistic amount of variability in the sampling distribution. The following R code implements this simulation, given K and ’level’:

r <- mu/(1 - mu)
B <- K/(1+r)
A <- r*B
aa <- A + 1
bb <- B + 1
rb <- rbeta(20, shape1=aa, shape2=bb)
