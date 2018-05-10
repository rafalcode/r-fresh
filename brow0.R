#!/usr/bin/env Rscript
fcn <- function(x, y){
  z <- x*y
  browser()
  return(z*exp(z))
}

val <- fcn(2, 1.3)
