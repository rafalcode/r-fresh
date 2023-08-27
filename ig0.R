#!/usr/bin/env Rscript
# this script does what? Using igraph!
library(igraph)
library(Cairo)

g <- make_ring(10)

set_vertex_attr(g, "name", value = letters[1:10])

# Cairo image template
CairoPNG("igplot0.png", 800, 800)
plot(g)
dev.off()

gz <- make_graph('Zachary')

gn <- make_graph(~ Alice-Bob:Claire:Frank, Claire-Alice:Dennis:Frank:Esther,
                George-Dennis:Frank, Dennis-Esther)
CairoPNG("igplot2.png", 800, 800)
plot(gn)
dev.off()

# So what's unusually about igraph is that the function syntax
# with V() and E() actually sets attributes, and need not be stored 
V(gn)$age <- c(25, 31, 18, 23, 47, 22, 50) 
V(gn)$gender <- c("f", "m", "f", "m", "m", "f", "m")
E(gn)$is_formal <- c(FALSE, FALSE, TRUE, TRUE, TRUE, FALSE, TRUE, FALSE, FALSE)
