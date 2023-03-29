#!/usr/bin/env Rscript
# this script does what?
args <- commandArgs(trailingOnly = TRUE)
numargs <- length(args)
enumargs <- 1 # expected numargs
if(numargs != enumargs) {
    print("This script converts and Rmd file to R, but retaining the documentation as comments.")
    print("It requires one argument, the Rmd file.")
    stop("Stopping right there.")
}

library(knitr)
purl(args[1], documentation=2)
