#!/usr/bin/env Rscript
# this script does what?
library(ggplot2)
library(Cairo)
library(PairedData)

p <- rpaired.gld(n=30,r=0.5)
# # put plot command here
data(lambda.table)

CairoPNG("pai0.png", 800, 800)
# p <- rpaired.gld(n=30,d1=lambda.table[7,],d2=lambda.table[7,],r=0.5)
pp <- plot(p)
show(pp)
dev.off()
