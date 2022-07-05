#!/usr/bin/env Rscript
# from Dirk Rcpp gallery
# http://dirk.eddelbuettel.com/blog/2012/09/02/
# I was confused abou this idea of binomial matrix
# well having run and inspected it, I can be more explicit
# this is about generating a random binary matrix, i.e. a matrix with random 1's and 0's.
library(inline)
# library(compiler) # for the cmpfun closure wrapper: not nec is recent versions of R
library(rbenchmark)

scott <- function(N, K) {
    mm <- matrix(0, N, K)
    apply(mm, c(1, 2), function(x) sample(c(0, 1), 1))
}
scottComp <- cmpfun(scott) # closure wrapping of some sort. Maybe it precompiles also?

# main
n <- 5
k <- 10
ss <- scottComp(n,k)
