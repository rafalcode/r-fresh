#!/usr/bin/env Rscript
# an example from a guy called Richard Morey.
library(Cairo)

# Distribution of p values
# Richard D. Morey
# 05/10/2016
# Distribution of t
# statistics (one sample t test)
# 
# The distribution of t
# statistics under an alternative hypothesis with true effect size δ and sample size n is
# t∼noncentral t_n-1(δn−−√).

n <- 50
delta <- .3
sigma <- 1

# simulation
v <- rnorm(n, sigma*delta, sigma)
tt0 <- t.test(v) # if you down specify mu (the given mean) then it will be set to zero.
# (the ones-sample t-test definitely needs mu, but in above command ts defautl is taken to be zeroQQZZ
# 
tstats <- replicate(10000, tt0$statistic)

## theoretical
tt <- seq( min(tstats), max(tstats), len = 250 )
yy <- dt(tt, df = n - 1, ncp = delta * sqrt(n) )

## plot
## simulation
CairoPNG("mor.png", 800, 800)
hist(tstats, freq = FALSE, col="gray", breaks = 20)
## theoretical
lines(tt, yy, lwd=2,col = "red")
dev.off()
