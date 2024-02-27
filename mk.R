#!/usr/bin/env Rscript
# this script does what? methylKit sim example
library(ggplot2)
library(Cairo)
library(methylKit)

# The methylation in 10% of the sites are elevated by 25%.
mymeth <- dataSim(replicates=4,sites=2000,treatment=c(1,1,0,0),
                   percentage=10,effect=25)
# the .Data componenet is a list of 16, 12 of which are 2000 integers.
# which vary from 0 to 700 though median is 24
# dunno.
stop("o")
CairoPNG("mk0.png", 800, 800)
gp <- ggplot(df, aes(x=x1, y=x2)) +
    geom_point(alpha=.2, size=3)
show(gp)
dev.off()

# scatter plots colorimetrics: I already hit up on it, they end up being symmetrical
