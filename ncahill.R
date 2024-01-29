#!/usr/bin/env Rscript
# this script does what?
# from one of rpubs thingies: https://rpubs.com/ncahill_stat/889056
# Niamh Cahill the author 
# Cannot make head or tail of it.
#
# Well, hang on. A very key point is that there is no independence between the various points.
# so you need function to define the nature of that independence
# for time points and considering a quantity such as hunger or some other wasting quantity
# the exp(-etc) function is good.

# although iid is hard to get, it's easy .. no functions or anything

library(R2jags)
library(runjags)
library(tidyverse)
library(tidybayes)
library(fields)
library(mvtnorm)
library(Cairo)

set.seed(28061989)

ntmpts <- 100 # number of time points.

# So, each time point is a separate random variable

yrspan <- 200 # span of years.
year <- sample(1:yrspan, ntmpts) %>% sort # unevenly spaced but ordered.
x <- year/ntmpts

# She talks about an autocorrelation function: it looks alot like a covariance matrix.
# the topic is Gaussian Processes, but her interest is time points, UNevenly spaced.

# Step 1: We’ll create a matrix of distances between every combination of time points using the rdist function from the fields package.

d <- fields::rdist(x) # a square matrix ntmpts x ntmpts
# actually shows the different between each time point

# Step 2: We’re going to create an autocorrelation function which we’ll call K. This will be a function of the distances in dist and we’re going to make the correlation function decay exponentially.
# NOTE that if they were EVENly spaced a function may not be necessary
# because the distance between the time points is fixed.
# however here that is not the case. SO THAT's WHY you need the function.
K <- exp(-d^2)
# looks alot like a covariance matrix.
# and actually that is what it's becoming. Basically we have a quantity that varies over time
# 

# To illustrate this, we’ll consider the first time point and we’ll look at how the correlation function decays as the distance between the time points increases.
# 
# plot(d[,1], K[,1], 
#      type = "l",
#      xlab = "year difference (hundred years)",
#      ylab = "autocorrelation")

# let's speed up decay ...
phi <- 5
K <- exp(-(phi^2)*(d^2))
# plot(d[,1], K[,1], 
#           type = "l",
#                xlab = "year difference (hundred years)",
#                ylab = "autocorrelation")

g <- mvtnorm::rmvnorm(5, sigma=K) # mean defaults to 0
# so this actually is 5 variates each
# matplot is part of standard graphics package
# matplot(x, t(g), 
#         type="l", 
#         ylab="g")
# that's actually quite a cute graph  but actually
# there needn't be similarity in the shapes at all
# but actually the key thing is the smoothness of the shapes.
# usually rnorms would be very spikey when you took lots of variates of each one.
# here we're only taking five but the important points is how each variate from 1 single drawing
# are dependant with each other.

# amplify the variation with:
sigma_g <- 2
g <- rmvnorm(5, sigma=(sigma_g^2)*K)

# BUT, actually what we'll do is stick to one variate for each variable/tmiepoint
# and add noise
sigma_y <- 0.3
eps <- rnorm(ntmpts, 0, sigma_y) # random noise

g <- rmvnorm(1, sigma=(sigma_g^2)*K) # GP

# additive noise
y <- c(g) + eps # Simulated data = GP + random noise

# for plotting
# BEWARE this is not curving fitting. If it were, it would be great.
# and looks fab.
# BUT it's not: noise is add to the "curve".
dat = tibble(year = year,
             x = x,
             y = y, 
             gp = c(g))

ggplot(dat, aes(x = year, y = y)) +
    geom_point(aes(colour = "Simulated Data")) +
    geom_line(aes(x = year, y = gp, colour = "Gaussian Process")) +
    labs(colour = "")

# She does go on howvever to fit .. but using JAGS
# In nay case the smoothed curve is what she calls a GP, and it's pretty
# much a flip of what I thought, what I thought were variates are actually separate random variables
