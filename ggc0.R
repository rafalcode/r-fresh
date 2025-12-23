#!/usr/bin/env Rscript
# this script does what? ggcX a seiries of scripts testin gout how ggplot used colours.
# the key here is 
# scale_fill_identity()
# it's the way of controlling colours that otherwise ggplot will do as it wants.
library(ggplot2)
library(Cairo)

# my colour plan - 17 colours.
mycols <-c("#A6CEE3", "#1F78B4", "#B2DF8A", "#33A02C", "#FB9A99", "#E31A1C", "#FDBF6F", "#FF7F00", "#CAB2D6", "#6A3D9A", "#FFFF99", "#B15928", "#F0027F", "#666666", "#FFFF33", "#8DD3C7", "#E6AB02")
stroke <-c("#333399", "#009900", "#FF6633", "#CCCCCC")

col0 <- c("red", "green", "blue", "yellow")

n <- 4 # how many different colours?
col0 <- mycols[1:n]
col0 <- stroke

df <- data.frame(
  x = 1:n,
  y = 1:n,
  Colour = col0)


CairoPNG("ggc0.png", 800, 800)
ggp <- ggplot(df, aes(x, y)) + geom_tile(aes(fill = Colour)) +
       scale_fill_identity()
show(ggp)
dev.off()
