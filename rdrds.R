#!/usr/bin/env Rscript
# this script does what? Read in an rds file and output as csv
args <- commandArgs(trailingOnly = TRUE)
numargs <- length(args)
enumargs <- 1 # expected numargs
if(numargs != enumargs) {
    print("This script retrieves and processes an R-generated RDS file")
    print("It requires one argument, the .rds file name (in full)")
    stop("Stopping right here")
}

# note rds coudl be anywhere (i.e. will have filepath,
# but csv is output in CWD (i.e. basename() is used)
fout <- gsub("(.+)\\.rds$", "\\1.csv", basename(args[1]))
rdscontent <- readRDS(args[1])
write.csv(rdscontent, fout, quote=F)
