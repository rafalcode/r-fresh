#!/usr/bin/env Rscript
# this script does what?
args <- commandArgs(trailingOnly = TRUE)
numargs <- length(args)
enumargs <- 2 # expected numargs
if(numargs != enumargs) {
    print("This script chunkifies an integer into fixed chunk sizes, with a remainder")
    print("1st arg should be the integer, 2nd arg is the chunk size")
    # warning("Stopping right here")
    stop("Stopping right here")
}
a1 <- as.integer(args[1])
a2 <- as.integer(args[2])

nfchunks <- floor(a1/a2) # number of full chunks.
remchunk <- a1%%a2 # the size of the remainder chunk

if(nfchunks) {
    for(i in 1:nfchunks) {
        for(j in 1:a2) {
            jj <- (i-1)*a2+j
            cat("chunk:", i, "/value:",jj,"\n")
        }
    }
}

if(remchunk) {
    for(j in 1:remchunk) {
        jj <- nfchunks*a2+j
        cat("chunk:", nfchunks+1, "/value:",jj,"\n")
    }
}
