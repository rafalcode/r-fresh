#!/usr/bin/env Rscript
# randomly mess with Std Dev of some genes
library(ggplot2)
library(Cairo)


ng <-20
bx <-4 # base expression.

df <- data.frame(BaseExp=rep(bx,ng), SD1=rnorm(ng,bx,sd=1), SD2=rnorm(ng, bx,sd=2), SD3=rnorm(ng,bx,sd=4))
df$A1 <- df$BaseExp+rnorm(ng,0,sd=1)
df$A2 <- df$A1 + rnorm(ng,0,sd=1)
df$A3 <- df$A2 +rnorm(ng,0,sd=1)
df$A4 <- df$A3 +rnorm(ng,0,sd=1)

#What was this?
# We see how progressive additions of sd=1 variates (which obviously is a liek a random walk)
# are ntohing compared to an sd=4, where there is wild variation.
# a power law at work witht he sd=4's!
