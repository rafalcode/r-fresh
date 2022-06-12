#!/usr/bin/env Rscript
# from:
# https://www.biostars.org/p/359307/#405010

library(Cairo)
library(ggplot2)

data <- read.csv("S1.csv", sep =";", header = TRUE, stringsAsFactors = FALSE)

CairoPNG("leite.png", 800, 800)
gg <- ggplot(data, aes(x= Group, y=Pathways, size=DEGs, color=FDR, group=Group)) +
    geom_point(alpha = 0.8) + 
    theme_classic() +
    scale_color_gradient(low = "red2",  high = "mediumblue", space = "Lab", limit = c(0.000000000000000007, 0.002)) +
    scale_size(range = c(2, 8))
show(gg)
dev.off()
