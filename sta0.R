#!/usr/bin/env Rscript

# how to get a static variable in R?
# see 

make.f <- function() {
    count <- 0
    f <- function(x) {
        count <<- count + 1
        return( list(mean=mean(x), count=count) )
    }
    return( f )
}

f1 <- make.f()
result <- f1(1:8)
print(result$count, result$mean)
result <- f1(1:8)
print(result$count, result$mean)

f2 <- make.f()
result <- f2(1:8)
print(result$count, result$mean)
result <- f2(1:8)
print(result$count, result$mean)
