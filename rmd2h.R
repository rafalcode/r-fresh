#!/usr/bin/env Rscript
# this script does what?
args <- commandArgs(trailingOnly = TRUE)
numargs <- length(args)
enumargs <- 1 # expected numargs
if(numargs != enumargs) {
    print("This script converts an Rmd fie to a html via knitr")
    print("It requires one argument")
    warning("Stopping right here")
}

library(knitr)
library(rmarkdown)

# knit2html(args[1])
rmarkdown::render(args[1])
