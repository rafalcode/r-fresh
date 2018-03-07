#!/usr/bin/env Rscript
# here we use an argument to form a new directory name and create it, then even outp a vector into as a file into that directory.
args <- commandArgs(trailingOnly = TRUE)
numargs <- length(args)
enumargs <- 1 # expected numargs
if(numargs != enumargs) {
    ostr<-sprintf("Usage: %i argument(s) required (hint: a distinct number for the creation of an output directory)", enumargs)
    print(ostr)
    stop("Bailing out ..")
}
oDir<-paste("nDir", toString(args[1]), sep="")
dir.create(oDir)
outF<-"outFle"
a<-seq(2,20,3)
write.csv(a, paste(oDir, "/", outF, sep=""))
