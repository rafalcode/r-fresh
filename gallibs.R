#!/usr/bin/env Rscript
# this script does what?
# gallibs.R, get all libs, to be run by sudo
# args <- commandArgs(trailingOnly = TRUE)
# numargs <- length(args)
# enumargs <- 1 # expected numargs
# if(numargs != enumargs) {
#     print("This script retrieves and processes a Gene Expression Omnibus entry")
#     print("It requires one argument, the accession serial number, eg: GSE49577")
#     warning("Stopping right here")
# }

args <- c("/opt/R-4.3.1")
# cat(paste0(args[1], "/bin\n"))
# cat(paste0(args[1], "/lib/R/library\n"))
libd <- paste0(args[1], "/lib/R/library")
# ld <- list.dirs(path=libd, recursive=F)
ld <- list.dirs(path=libd, recursive=F, full.names=F)
# ld <- gsub("
# cat(paste(ld, collapse="\n"))
BiocManager::install(ld)
