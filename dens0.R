#!/usr/bin/env Rscript
# this script does what? let's have a look at density func ... it used by a lot of things (oh say svaseq?)
library(ggplot2)
library(Cairo)

CairoPNG("dens00.png", 800, 800)
plot(density(c(-20, rep(0,98), 20)), xlim = c(-4, 4))  # IQR = 0
dev.off()

# The Old Faithful geyser data
d <- density(faithful$eruptions, bw = "sj")
CairoPNG("dens01.png", 800, 800)
plot(d)
dev.off()

CairoPNG("dens02.png", 800, 800)
plot(d, type = "n")
polygon(d, col = "wheat")
dev.off()

## Missing values:
x <- xx <- faithful$eruptions
x[i.out <- sample(length(x), 10)] <- NA
doR <- density(x, bw = 0.15, na.rm = TRUE)
lines(doR, col = "blue")
points(xx[i.out], rep(0.01, 10))

## Weighted observations:
fe <- sort(faithful$eruptions) # has quite a few non-unique values
## use 'counts / n' as weights:
dw <- density(unique(fe), weights = table(fe)/length(fe), bw = d$bw)
utils::str(dw) ## smaller n: only 126, but identical estimate:
stopifnot(all.equal(d[1:3], dw[1:3]))

## simulation from a density() fit:
# a kernel density fit is an equally-weighted mixture.
fit <- density(xx)
N <- 1e6
x.new <- rnorm(N, sample(xx, size = N, replace = TRUE), fit$bw)
plot(fit)
lines(density(x.new), col = "blue")

# Cairo image template
CairoPNG("dens00.png", 800, 800)
# put plot command here
dev.off()
