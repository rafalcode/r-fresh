#!/usr/bin/env Rscript
# this script does what?
library(ggplot2)
library(Cairo)
library(qvalue)

data(hedenfalk)

pvalues <- hedenfalk$p
qobj <- qvalue(p = pvalues)



# Cairo image template
# CairoPNG("fname.png", 800, 800)
# put plot command here
# dev.off()
