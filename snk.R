#!/usr/bin/env Rscript
# this script does what?
library(ggplot2)
library(Cairo)
library(networkD3)

nodes = data.frame("name" = c("Node A", "Node B", "Node C", "Node D"))

links = as.data.frame(matrix(c( 0, 1, 10, # Each row represents a link. The first number
                                0, 2, 20, # represents the node being conntected from. 
                                1, 3, 30, # the second number represents the node connected to.
                                2, 3, 40),# The third number is the value of the node
                              byrow = TRUE, ncol = 3))

names(links) = c("source", "target", "value")

sankeyNetwork(Links = links, Nodes = nodes,
               Source = "source", Target = "target",
                Value = "value", NodeID = "name",
                fontSize= 12, nodeWidth = 30)

# Cairo image template
# put plot command here
# CairoPNG("fname.png", 800, 800)
# dev.off()
