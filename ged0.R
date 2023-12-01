#!/usr/bin/env Rscript
# this script does what?
# another attempt to sort out grid grobs and gedit
library(ggplot2)
library(Cairo)
library(grid)

my_data <- data.frame(x = sample(3592185:3850502, 10), y = sample(4104561:4258891, 10), var = c("A","A","B","A","A","A","B","A","B","A"))

g <- ggplot() + 
    geom_point(data=my_data, aes(x = x, y = y, colour = var, shape = var, fill = var), stroke = 0.4) +
    scale_shape_manual("",values=c(24,24)) +
    scale_color_manual("",values=c("gray30", "steelblue")) +
    scale_fill_manual("",values=c("gray90", "lightblue")) +
    geom_point(data=my_data, aes(x = x, y = y, fill = var), color = "black",  shape = 16, size = 0.2) +   
    guides(fill = guide_legend(override.aes = list(size = c(8,8)))) +
    labs(title = NULL, x = NULL, y = NULL)  

# CairoPNG("ged0.png", 800, 800)
# show(g)
# dev.off()

grid.ls(grid.force())    # get the names of all the grobs in the ggplot

# The edit - to set the size of the point in the legend to 1.5 mm
#grid.gedit("key-[-0-9]-1-1", size = unit(1.5, "mm"))    # 1-1 interacts with the symbols from the first geom point (triangles)
grid.gedit("key-[-0-9]-1-2", size = unit(1.5, "mm"))    # 1-2 interacts with the symbols from the second geom point (dots)

gg <- grid.grab()
ggsave(plot=gg, file="test_plot.png", units="in", width=8.8, height=6, dpi = 300, type ="cairo-png")
