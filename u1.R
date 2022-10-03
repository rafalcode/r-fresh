#!/usr/bin/env Rscript
# this script does what?
# bit of reverse eng. as of ocurse no data file given, but some output yes.
# ref. https://rstudio-pubs-static.s3.amazonaws.com/443636_b649615b525f437bad4ac9fa53ebe56d.html
# only worth it for the start

# no data, so simulate 50 rows of 2 vars each
tab <- read.table("u1.tab")
m1 <- mean(tab$V1)
m2 <- mean(tab$V2)
sd1 <- sqrt(var(tab$V1))
sd2 <- sqrt(var(tab$V2))

nr <- 50 # one supposes.
hh <- 2 # height for cutree

x <- data.frame(X1=rnorm(nr, mean=m1, sd=sd1), X2=rnorm(nr, mean=m2, sd=sd2))

# It should be clear, the pairwise distance matrix is the key
dx <- dist(x)
hclust.out <- hclust(dx)

# Cut by height
# cu <- cutree(hclust.out, h = hh) # cut tree by height
cu <- cutree(hclust.out, k=4) # by number of clusts.

plot(hclust.out)
abline(h = hh, col = "red")
