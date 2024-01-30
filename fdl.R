#!/usr/bin/env Rscript
# this script does what?
# actually trying for quite a useful utility
# sho repeated whole lines in a file
# args <- commandArgs(trailingOnly = TRUE)
# numargs <- length(args)
# enumargs <- 1 # expected numargs
# if(numargs != enumargs) {
#     print("fdl.R find duplicate lines in a file, 1 arg: the file\n")
#     stop("Stopping right here")
# }

# rl <- readLines(args[1])
# rl <- readLines("sim_roc_rf0.R")
rl <- readLines("t.a")

uni <- unique(rl)

w <- match(rl, uni)
wt <- table(w)
ww <- unname(which(wt>1))
wt2 <- sort(wt[ww], decreasing=T)
for(i in 1:length(wt2)) {
        rli <- as.integer(names(wt2)[i])
        rtimes <- unname(wt2[i])
        cat(paste0("Times ", rtimes, ": ", rl[rli], "\n"))
}
