#!/usr/bin/env Rscript
# this script does what?
# this webpage https://ahoyyangbai.wordpress.com/2013/08/23/simple-brownian-motion-generation-in-r/
# came up in relation to Wiener processes a  generalization of a brownian motion
# which is just a simple random walk.
# that webpage is quite bad though

library(Cairo)
set.seed(1)

# 1. rn<-rnorm(n=1000,m=0,sd=1)
# well, m=0, sd=1 is the default so zap that.
rn <- rnorm(1000)

CairoPNG("simpbrow0.png", 800, 800)
plot(rn, type='l')
dev.off()

CairoPNG("simpbrow1.png", 800, 800)
hist(rn,breaks=25)
dev.off()

# so a random walk aka brownian mot ... simply add one onto the other .. that simple
# what's dp ... dipslacment? Doubt it. Delt postion? They're not deltas.
# quality of this shit article doesn't merit decipyhering this.
dp <- cumsum(rn)

CairoPNG("simpbrow2.png", 800, 800)
plot(dp,type='l')
dev.off()

# And then this crazy step size digression
# which is senseless ... change step size ... a silly endeavour, does not advance learning
# If anyone want to make the x-axis show a certain step size, just assign values to two vectors and plot one against the other. Codes are like:

# another subopt statement
# x<-c(1:100)
x <- 1:100

x<-x/100##adjust the step size of my 100 observations to be 0.01. Sheesh!

y<-rnorm(100, 0, 0.01)##100 observations here with variance 0.01

y<-cumsum(y)

CairoPNG("simpbrow3.png", 800, 800)
plot(x, y, type='l', xlab='Time', ylab='Value')
dev.off()

# what a pointless exercise.
