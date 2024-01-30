#!/usr/bin/env Rscript
# this script does what?
# a walk through
# ref. https://jeremydfoote.com/Communication-and-Social-Networks/resources/ggraph_walkthrough.html
library(Cairo)
library(igraph)
library(tidygraph)
library(ggraph)
library(tidyverse)

set_graph_style() # This sets the default style to the graph style

G0 <- erdos.renyi.game(50, .4)
G <- as_tbl_graph(G0)

