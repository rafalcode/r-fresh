#!/usr/bin/env Rscript
# this script does what? checks arguments
args <- commandArgs(trailingOnly = TRUE)
numargs <- length(args)
enumargs <- 1 # expected numargs
if(numargs != enumargs) {
    print("This script retrieves and processes a Gene Expression Omnibus entry")
    print("It requires one argument, the accession serial number, eg: GSE49577")
    warning("Stopping right here")
}
str(args[1])
