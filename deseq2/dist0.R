#!/usr/bin/env Rscript
# this script does what?
# from: https://stackoverflow.com/questions/9879608/how-do-i-manipulate-access-elements-of-an-instance-of-dist-class-using-core-r
library(ggplot2)
library(Cairo)

# Cairo image template
# CairoPNG("fname.png", 800, 800)
# put plot command here
# dev.off()

distidx <- function(i,j,n)
{
    #given row, column, and n, return index
    n*(i-1) - i*(i-1)/2 + j-i
}

rowcol<-function(ix,n) #given index, return row and column
{
    # less useful, this guy, in my opinion.
    # given the row and column, tells you which points it refers to
    # alot less use for this.
    nr=ceiling(n-(1+sqrt(1+4*(n^2-n-2*ix)))/2)
    nc=n-(2*n-nr+1)*nr/2+ix+nr
    cbind(nr,nc)
}

# A little test harness to show it works:
set.seed(42)
npts <- 8
testd <- dist(rnorm(npts))
# that'a distance object a lower triangular
# it's easy to turn it into a fully fledged matrix with diag as zeros:
dasm <- as.matrix(testd)

# So you can ask a certain distance between two points from it
# say the second from the fourth
# dasm[2,4]   #row<col
# you'll get sme value if you reverse the numbers
# because that matrix carries double the redundancy

# but say you don't want that, you want to work witht he dist object
# so what is the index where the distance between 2 and 4 can be found?
# using distidx() function you have: 
distidx(2, 4, npts) 
# so
testd[distidx(2, 4, npts)]
# will give you the value.

# So generally we can see that single indexing works in usual R fashion, Column Major.

# testd[c(42,119)]
rowcol(c(2,7), npts)  # = (3,8) and (8,15)
as.matrix(testd)[3,8]
# as.matrix(testd)[8,15]
