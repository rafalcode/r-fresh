#!/usr/bin/env Rscript
# this script does what? checks arguments
args <- commandArgs(trailingOnly = TRUE)
numargs <- length(args)
enumargs <- 0 # expected numargs
if(numargs != enumargs) {
    print("This script takes no arguments")
    warning("Stopping right here")
    # won't exit though.
}
va <- c(1095, 312, 1950, 1044, 1749, 621, 1650, 2067, 999, 1437, 753, 1263, 73, 2118, 1668, 2709, 4860, 1725)
print(length(va))
