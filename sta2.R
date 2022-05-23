#!/usr/bin/env Rscript

# For static variables as used in C, R
# needs to resot to closures, which ar efunct

loopidx <- function() {
  a <- runif(4)
  i <- 0
  sz <- length(a)
  function() {
    # do something useful, then ...
    i <<- i%%sz + 1 # a local-to-parent assignment
    return(a[i])
  }
}

me <- loopidx()
for(i in 1:7) {
    # cat(paste0("No. ", i, " is ", me(), "\n"))
    print(me())
}
