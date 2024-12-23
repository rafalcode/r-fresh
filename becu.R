#!/usr/bin/env Rscript
# generate beta distribution curves
# NOTE this is not he same 
library(ggplot2)
library(Cairo)

args <- commandArgs(trailingOnly = TRUE)
numargs <- length(args)
enumargs <- 2 # expected numargs
if(numargs != enumargs) {
    print("This script generates a beta distribution density curve, needing two arguments, shape1 and shape2")
    stop("Stopping right here")
}

sh1 <- as.numeric(args[1])
sh2 <- as.numeric(args[2])

sq <- seq(0, 1, by=0.005)
bvals <- dbeta(sq, sh1, sh2) # uniform beta vals (the x is uniform, not the y!)
# above comment is impotant, because when we get narrow distri shapes, the dot representation is not visually satisfying
# because very few dots gather aroudn the narrow part!

CairoPNG("db0.png", 800, 800)
# lines connecting the dots (the type value is key here)
# plot(sq, bvals, pch=20, cex=1.5, col="firebrick", type="b")
plot(sq, bvals, pch=20, cex=1.5, col="firebrick")
# lines(sq, bvals, col="grey")
# with lines you have more control over colour of the lines
dev.off()


# Ques, how about instead of dots we have density type line graph?
# not lines(), thsat just connects the dots. A smoothing curve I suppose.
# A density on  lots of rand variates would work of course.
# any alternatives?

# Actually this is a common issue, sort of surprising I never hit it before,
# at least in this guise.
# so definitely, this is about smoothing.
# say you only have 20 points from your dbeta(),. For this exercise
# things are very clear, the smoothing will work
# and we don't particularly need it, but let's pretend it's real data
# and we want something curvy coming out.
# of course we're entering into splines and stuff here.
# but this seems a good post:
# https://stackoverflow.com/questions/3480388/how-to-fit-a-smooth-curve-to-my-data-in-r
# "loess()" is veryone's favourite", look no further!

# incredibly enough, I discounted lines() above, becuase I thought there would be somehting more sophiticated happening.

# however, the winner, is not to use lines, but the type="l" which is nativr to plot.
# You'll need probably more than 100 points, 1000 say, but that's still good
# and you get a very neat curve.
lo2 <- loess(bvals~sq)
print(str(lo2))
CairoPNG("db1.png", 800, 800)
plot(sq, bvals, pch=20, cex=0.5, col="steelblue", type="l")
# lines(predict(lo2), col='black', lwd=2)
dev.off()


# but this does work for Dirk:
# but actually it's not spliend between points.
x <- 1:10
y <- c(2,4,6,8,7,12,14,16,18,20)
lo <- loess(y~x)
print(str(lo))

CairoPNG("dirk0.png", 800, 800)
plot(x,y)
lines(predict(lo), col='steelblue', lwd=2)
dev.off()
