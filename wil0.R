#!/usr/bin/env Rscript
# The wikipedia example: it get p-value 0.6113
# having ties and zero really complicates it .. not a very good example!
library(ggplot2)
library(Cairo)

raw <- c(125, 110, 115, 122, 130, 125, 140, 120, 140, 140, 115, 124, 140, 123, 125, 137, 140, 135, 135, 145)
lr <- length(raw)
x1 <- raw[seq(1,lr,2)]
x2 <- raw[seq(2,lr,2)]

wilc <- wilcox.test(x1, x2) # I don't get the same, something close.
# the R function doesn't like zero and ties
# Cairo image template
# CairoPNG("fname.png", 800, 800)
# put plot command here
# dev.off()
