#!/usr/bin/env Rscript
# this script does what?
# working but it's not doing what I want.
library(igraph)
library(ggraph)
library(ggplot2)
library(Cairo)

col_fun <- colorRampPalette(c('tomato', 'skyblue'))

set.seed(1)
g <- erdos.renyi.game(10, .05)
V(g)$name <- letters[1:10]
V(g)$cap <- 1:10
idx <- 1:10
# idx[c(92, 2)] <- c(2, 92)
# g <- permute(g, idx)
# V(g)$size <- scales::rescale(degree(g), c(3, 20))
# V(g)$color <- col_fun(vcount(g))
# V(g)$color[2] <- "#FF0000"
# V(g)$color[9] <- "#00FF00"
# V(g)$color <- ifelse(V(g)$name<6, "yellow", "blue")
V(g)$color <- ifelse(V(g)$cap<6, "tomato", "skyblue")
V(g)$size <- ifelse((V(g)$cap%%2)==0, 24, 48)

# Cairo image template
p <- ggraph(g, layout="fr")
p <- p + ggraph::geom_node_point(ggplot2::aes(fill=.data$color, size=.data$size), shape=21, color="black")
p <- p + geom_node_text(aes(label = name), repel=TRUE)
CairoPNG("erdos.png", 800, 800)
plot(p)
dev.off()
