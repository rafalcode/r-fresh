#!/usr/bin/env Rscript
# I want to find the sequence that has exactly 2 < -2 and 2> 2

for(i in 1:100) {
    rx <- sort(rnorm(100))
    if((rx[1]<(-2)) & rx[2]<(-2) & (rx[3]>(-2)) & (rx[98]<(2)) & (rx[99]>(2)) & (rx[100]>(2))) {
        # cat(paste0("Yes for i=", i), "\n")
        cat("i=", i, " - ", paste(rx[c(1:3,98:100)], collapse=" "), "\n")
    }
}
