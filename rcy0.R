library(RCy3)
library(igraph)

ig <- make_graph("Zachary")
createNetworkFromIgraph(ig)

# you get 34 vertices and 76 or so edges.
# there's no names on this one (no Mr. H)
# you just get the vertex indices.
