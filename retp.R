#!/usr/bin/env Rscript
# this script does what?

rd9 <- function(n, frac) # dodgy random sequence: enter total number and fraction to be 9's.
{
    s <- as.integer(runif(n)*10)
    s[sample(1:n, round(frac*n))] <- 9
    s
}

cat(paste0("0.", paste0(rd9(16,0.5), collapse="")), "\n")
