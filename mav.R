#!/usr/bin/env Rscript
# ref. Neil 2015 question 
library(microbenchmark)
set.seed(1)
ma0 <- function (x , k ) {
   n <- length ( x ) # length of input vector
   nMA <- n - k +1 # the number of elements in the moving average
   xnew <- rep ( NA , nMA ) # to store moving average
   i <- 1 # counter variable
   # to calculate moving average , will calculate average of
   # the next k elements starting from position i
   while ( i <= nMA ){ # until the last moving average
      # calculate average for k values starting from element i
      xnew [ i ] <- mean ( x [ ( i :( i +k -1) ) ])
      i <- i +1
   }
   xnew
}

ma1 <- function (x , k ) {
   n <- length ( x ) # length of input vector
   nMA <- n - k +1 # the number of elements in the moving average
   xnew <- sapply(1:nMA, function(y){sum(x[y:(y+k-1)])/k})
   xnew
}

ma2 <- function (x , k ) {
   n <- length ( x ) # length of input vector
   xbk <- x/k
   nMA <- n - k +1 # the number of elements in the moving average
   xnew <- sapply(1:nMA, function(y){sum(xbk[y:(y+k-1)])})
   xnew
}
 
# x <- c (3.5 , 3.2 , 2.9 , 3.1 , 2.9 , 2.8 , 3.0 , 2.7)
x <- runif(10000)

microbenchmark (
ma1(x, k =5),
ma2(x, k =5)
)
