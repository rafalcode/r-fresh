#!/usr/bin/env Rscript
# this script does what?
# martingales have a great name, why don't I know more about them
# check this one:
# https://stackoverflow.com/questions/69419602/how-to-simulate-a-martingale-process-problem-in-r
library(Cairo)
args <- commandArgs(trailingOnly = TRUE)
numargs <- length(args)
enumargs <- 2 # expected numargs
if(numargs != enumargs) {
    print("This script imulates a suppsed martingale\n")
    print("It requires 2 arguments: 1) seed 2) quantity of cloakrooms guests\n")
    stop("Stopping right here")
}

set.seed(args[1])
n <- args[2]

a <- 1:n
s <- sample(a, n)

cat("First sequence is indexed coat owners, second is random allocation, third are the indxed that are correct and will be deleted in next iter\n")
# for (i in 1:6) {
times <- 0
while(length(a)>0) {
    cat(paste0(a, collapse=" "))
    cat("\n")
    cat(paste0(s, collapse=" "))
    cat("\n");
    w <- which(a == s)
    if(length(w)>0) {
        a <- a[-w]
        cat(paste0(w, collapse=" "))
    } else {
        cat("No person gets right coat!")
    }
    cat("\n");
    cat("\n");
    s <- sample(a, length(a))
    times <- times +1
}
cat(paste0("Took ", times, " iterations.\n"))
