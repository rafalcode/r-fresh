#!/usr/bin/env Rscript
## this script does what?
library(Cairo)
library(RUVcorr)

set.seed(400)
Yind <- simulateGEdata(n=3000, m=1000, k=10, size.alpha=2,
                       corr.strength=5, g=NULL, Sigma.eps=0.1, 
                       nc=2000, ne=1000, intercept=TRUE, check=TRUE)
# n, nunber of genes
# m, num samples.





# Cairo image template
# CairoPNG("fname.png", 800, 800)
# put plot command here
# dev.off()
