#!/usr/bin/env Rscript
# this script does what?
library(ggplot2)
library(Cairo)


# Illutration of Berksons paradox

bq <- .2 # beauty quotient, 1 or 5 people can be rated "attractive"
sq <- .1 # singing quotient, 1 in 10 are good singers.

n <- 10000 # size of pool of people

df <- data.frame(apv=ifelse(runif(n)>bq, "notatt", "yesatt"),
                 spv=ifelse(runif(n)>sq, "badsing", "goodsing"))

wa <- which(df$apv=="yesatt")
