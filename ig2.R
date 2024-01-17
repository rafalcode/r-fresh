#!/usr/bin/env Rscript
# this script does what? From the drawring raphs part of the manual.
library(igraph)
# args <- commandArgs(trailingOnly = TRUE)
# numargs <- length(args)
# enumargs <- 1 # expected numargs
# if(numargs != enumargs) {
#     print("This script retrieves and processes a Gene Expression Omnibus entry")
#     print("It requires one argument, the accession serial number, eg: GSE49577")
#     warning("Stopping right here")
# }

# plotting a simple ring graph, all default parameters, except the layout
g <- make_ring(10)
g$layout <- layout_in_circle

# three rendering engines:
# plot(g)
# tkplot(g)
# rglplot(g)

# plotting a random graph, set the parameters in the command arguments
g1 <- barabasi.game(100)
# plot(g, layout=layout_with_fr, vertex.size=4, vertex.label.dist=0.5, vertex.color="red", edge.arrow.size=0.5)

# plot a random graph, different color for each component
g2 <- sample_gnp(100, 1/100)
comps <- components(g)$membership
colbar <- rainbow(max(comps)+1)
V(g)$color <- colbar[comps+1]
# plot(g, layout=layout_with_fr, vertex.size=5, vertex.label=NA)

# plot communities in a graph
g3 <- make_full_graph(5) %du% make_full_graph(5) %du% make_full_graph(5)
g33 <- add_edges(g3, c(1,6, 1,11, 6,11))
com <- cluster_spinglass(g33, spins=5)
V(g33)$color <- com$membership+1
g33 <- set_graph_attr(g33, "layout", layout_with_kk(g33))
plot(g33, vertex.label.dist=1.5)

# draw a bunch of trees, fix layout
# igraph_options(plot.layout=layout_as_tree)
# plot(make_tree(20, 2))
# plot(make_tree(50, 3), vertex.size=3, vertex.label=NA)
# tkplot(make_tree(50, 2, mode="undirected"), vertex.size=10, vertex.color="green")

