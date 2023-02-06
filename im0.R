#!/usr/bin/env Rscript
# this script does what?
library(Cairo)

# CairoPNG("im0.png", 800, 800)
png(filename = "im0.png", width = 5, height = 15, units = "cm", pointsize = 12, res = 300)

par(bg = "transparent", usr = c(0, 51, 0, 451)) #make the plot window a certain size?
plot(x=NULL, y=NULL , type = "n", axes = F, xlab = "", ylab = "", xlim=c(0,51), ylim=c(0,450)) #set up the plot?

#draw rectangles for thermometer
rect(0, 0, 50, 148, col = "#c00000", border = "transparent")     #red
rect(0, 148, 50, 225, col = "#ed7d31", border = "transparent") #orange
rect(0, 225, 50, 297, col = "#ffc000", border = "transparent") #gold
rect(0, 297, 50, 360, col = "#92d050", border = "transparent") #lgreen
rect(0, 360, 50, 450, col = "#00b050", border = "transparent") #dgreen

dev.off()
