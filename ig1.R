# second script for ig1
# illustrating how to set nodes colours based on an attribute.
library(Cairo)
library(igraph)
library(RColorBrewer)

# create an example network
g <- make_ring(5)

# assign vertex attributes
g <- set.vertex.attribute(g, 'group', 1, 'A')
g <- set.vertex.attribute(g, 'group', 2, 'A')
g <- set.vertex.attribute(g, 'group', 3, 'B')
g <- set.vertex.attribute(g, 'group', 4, 'B')
g <- set.vertex.attribute(g, 'group', 5, 'C')

# create color pallet based on unique values for vertex attribute
pal <- brewer.pal(length(unique(V(g)$group)), "Dark2")

# plot network
CairoPNG("ringgra.png", 800, 800)
# plot(g, vertex.color = "gray")
plot(g, vertex.color = pal[as.numeric(as.factor(vertex_attr(g, "group")))])
dev.off()
