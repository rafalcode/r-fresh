# using the loess curve as covered here:
# https://stackoverflow.com/questions/3480388/how-to-fit-a-smooth-curve-to-my-data-in-r
# which I gather Dirk wrote.
library(ggplot2)
library(Cairo)

x <- 1:10
y <- c(2,4,6,8,7,12,14,16,18,20)
lo <- loess(y~x)

CairoPNG("loes0.png", 800, 800)
plot(x,y)
lines(predict(lo), col='red', lwd=2)
dev.off()

# the main thing to note here is that it does not go though all the points.
