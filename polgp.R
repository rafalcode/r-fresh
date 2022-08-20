#!/usr/bin/env Rscript
# A test of the polyester package's get_params for zero inf negbin distri's.
library(ggplot2)
library(Cairo)

library(matrixStats)
library(ballgown)
library(polyester)
library(RSkittleBrewer)
trop = RSkittleBrewer('tropical')

data(bg) # gives you bg, a ballgown object, fairly lengthy
countmat <- fpkm_to_counts(bg, mean_rps=400000)
params <- get_params(countmat) # get list of 4: p0, mu, size, fit.
lcounts <- log(countmat+1)

rm0 <- rowMeans(lcounts)
rv0 <- rowVars(lcounts)

CairoPNG("meanvarrel.png", 800, 800) # mean - variance relationship
# plot(rm0, rv0, pch=19,cex=0.25,col=trop[1],main="zebrafish data w/fit")
plot(rm0, rv0, pch=19, col=trop[1],main="zebrafish data w/fit")
lines(params$fit,col=trop[2])
dev.off()
# So this trick which tries to generate a line from the polyester params$fit
# is still odd ... I don't kno w what's its mapping.
