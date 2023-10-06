#!/usr/bin/env Rscript
# this script does what? I'm using chatgpt and it's very wrong.
# Load the packages
library(networkD3)
library(htmlwidgets)
library(webshot)
library(rsvg)
library(Cairo)

# Sample data
# Load data
data(MisLinks)
data(MisNodes)
# Create graph
sn0 <- forceNetwork(Links = MisLinks, Nodes = MisNodes, Source = "source",
             Target = "target", Value = "value", NodeID = "name",
             Group = "group", opacity = 0.4, zoom = TRUE)

# Create graph with legend and varying node radius
sn1 <- forceNetwork(Links = MisLinks, Nodes = MisNodes, Source = "source",
             Target = "target", Value = "value", NodeID = "name",
             Nodesize = "size",
             radiusCalculation = "Math.sqrt(d.nodesize)+6",
             Group = "group", opacity = 0.4, legend = TRUE)

# Create graph directed arrows
sn2 <- forceNetwork(Links = MisLinks, Nodes = MisNodes, Source = "source",
             Target = "target", Value = "value", NodeID = "name",
             Group = "group", opacity = 0.4, arrows = TRUE)

# Save the networkD3 graph to an HTML file
saveWidget(sn0, "sn0.html")
# or:
# webshot::save_svg(sn0, "sn0.svg")
# noe that function doesn't exist.
# rsvg::rsvg_pdf(sn0, "sn0.pdf")
# CairoSVG(sn0, "sn0.pdf")
CairoSVG(12, 12, "sn0.html", "sn0.pdf")
