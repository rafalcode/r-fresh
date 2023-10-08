#!/usr/bin/env Rscript
# this script does what?
# this is the rendering of 
# Neil Hurley's favourite:
# in matlab:
# N=10
# x=rand(N,1)
# y=rand(N,1)
# tri=delaunay(x,y)
# I asked chatgpt and it said use deldir.
library(ggplot2)
library(Cairo)
library(deldir)

set.seed(42)
x <- runif(20)
y <- runif(20)
d <- deldir(x,y)

ttt <- triang.list(d)

CairoPNG("delau0.png", 800, 800)
plot(ttt,border="red",showrect=TRUE,rectcol="green")
dev.off()

sss <- tile.list(d)

CairoPNG("delau1.png", 800, 800)
plot(sss)
plot(ttt, add=T, border="blue", showrect=T, rectcol="red") # that add superimposes on previous plot
dev.off()
