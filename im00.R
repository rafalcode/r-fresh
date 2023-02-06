#!/usr/bin/env Rscript
# this script does what?

png(filename = "im00.png", width = 5, height = 15, units = "cm", pointsize = 12, res = 300)

# par(bg = "transparent", usr = c(0, 51, 0, 451)) #make the plot window a certain size?
par(bg = "white", usr = c(0, 51, 0, 451)) #make the plot window a certain size?
plot(x=NULL, y=NULL , type = "n", axes = F, xlab = "", ylab = "", xlim=c(0,51), ylim=c(0,450)) #set up the plot?

#draw rectangles for thermometer
xx <- c(148, 225, 297, 360)
yy <- rep(50,4)

# points(148, 50, col="black", pch = 19, cex = 2)
# points(225, 50, col="black", pch = 19, cex = 2)
# points(297, 50, col="black", pch = 19, cex = 2)
# points(360, 50, col="black", pch = 19, cex = 2)
points(xx, yy, col="black", pch = 19, cex = 2)

dev.off()
