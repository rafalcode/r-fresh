#!/usr/bin/env Rscript

# For static variables as used in C, R
# needs to resot to closures, which ar efunct

loopidx <- function(arr) {
  i <- 0
  sz <- length(arr)
  function(arr) {
    # do something useful, then ...
    i <<- i%%sz + 1 # a local-to-parent assignment
    return(arr[i])
  }
}

a <- runif(4)
a
me <- loopidx(a)
for(i in 1:7) {
    res <- me(a)
    cat(paste0("No. ", i, " is ", res, "\n"))
}
