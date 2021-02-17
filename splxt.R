#!/usr/bin/env Rscript
# this script does what?
args <- commandArgs(trailingOnly = TRUE)
numargs <- length(args)
enumargs <- 1 # expected numargs
if(numargs != enumargs) {
    print("This script takes a file name and strips the etxension from it")
    stop("Stopping right here")
}

# split on dot
ss <- strsplit(args[1], "\\.")[[1]]
# but we may have more than one dot! If so the last one is the most important:
lss <- length(ss)
fss <- paste(ss[1:lss-1], collapse='.')
# OK fss holds eveythign until the last dot. So the new extension can now be:
newn <- paste0(fss, ".newextension")
newn
