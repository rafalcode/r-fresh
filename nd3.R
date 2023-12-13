#!/usr/bin/env Rscript
# this script does what?
# Load the packages
library(networkD3)
library(htmlwidgets)
library(webshot)
# from chatgpt actually, got tired of looking elsewhere.

# Sample data
nodes <- data.frame(name = c("A", "B", "C"))
links <- data.frame(source = c(0, 0, 1), target = c(1, 2, 2))

# Create a simple networkD3 graph
simpleNetwork <- forceNetwork(Links = links, Nodes = nodes, NodeID="name")

# Use the htmlwidgets package to view the networkD3 visualization in your R console or RStudio viewer pane:
# View the networkD3 graph
# simpleNetwork
# This will display the networkD3 visualization in your R console or RStudio viewer pane. You can interact with the graph and explore the nodes and links as needed.

# If you want to save the networkD3 visualization to an HTML file for sharing or embedding in a web page, you can use the saveWidget function from the htmlwidgets package:

# Save the networkD3 graph to an HTML file
saveWidget(simpleNetwork, file = "networkD3_graph.html")
# or:
webshot::save_svg(simpleNetwork, file = "nd3.svg")
