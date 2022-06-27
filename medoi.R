#!/usr/bin/env Rscript
# this script does what?
library(ggplot2)
library(Cairo)
library(cluster)

set.seed(1)
## generate 25 objects, divided into 2 clusters.
a <- rnorm(10,0,0.5)
b <- rnorm(10,0,0.5)
c <- rnorm(15,5,0.5)
d <- rnorm(15,5,0.5)
x <- rbind(cbind(a, b), cbind(c, d))

pamx <- pam(x, 2)
pamx # Medoids: '7' and '25' ...
# summary(pamx)
CairoPNG("pamx.png", 800, 800)
plot(pamx)
dev.off()

## use obs. 1 & 16 as starting medoids -- same result (typically)
(p2m <- pam(x, 2, medoids = c(1,16)))

## no _build_ *and* no _swap_ phase: just cluster all obs. around (1, 16):
p2.s <- pam(x, 2, medoids = c(1,16), do.swap = FALSE)
# p2.s
p3m <- pam(x, 3, trace = 2)

## rather stupid initial medoids:
(p3m. <- pam(x, 3, medoids = 3:1, trace = 1))
# pam(daisy(x, metric = "manhattan"), 2, diss = TRUE)
# data(ruspini)
## Plot similar to Figure 4 in Stryuf et al (1996)
## Not run: plot(pam(ruspini, 4), ask = TRUE)

