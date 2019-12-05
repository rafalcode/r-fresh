#!/usr/bin/env Rscript
# inovkes a system call to awk, using its -v option.

args <- commandArgs(trailingOnly = TRUE)
numargs <- length(args)
enumargs <- 1 # expected numargs
if(numargs != enumargs) {
    stop("This script takes 1 argument: a \".fam\" file.")
}


dateval <- "2019"
# sysCall <- paste0("awk -v var=\"", dateval, "\" '{$1=var\"_\"$1; print}' ", dset)
sysCall <- paste0("awk -v var=\"", dateval, "\" '{$1=var\"_\"$2; $2=var\"_\"$2; print}' ", args[1])

# cat(paste0(sysCall, "\n"))
rc <- system(sysCall, intern=FALSE)
if (rc != 0) {
    stop("sysCall didn't go well")
}
