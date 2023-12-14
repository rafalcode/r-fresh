library(ggraph)
library(igraph)
library(igraphdata)
library(tidyverse)
library(tidygraph) # for as_tbl_graph
# library(ggplot2)
# library(reshape2)
library(Cairo)

data(karate)
data(highschool)

graph <- as_tbl_graph(karate) %>% mutate(degree = centrality_degree())
# graph <- as_tbl_graph(highschool) %>% mutate(degree = centrality_degree())

lapply(c('stress', 'fr', 'lgl', 'graphopt'), function(layout) {
  gg <- ggraph(graph, layout = layout) + 
    # geom_edge_link this is for the edges only
    # geom_edge_link(aes(colour = factor(color)), show.legend = FALSE) +
    geom_edge_link(aes(colour = factor(weight)), show.legend = FALSE) +
    geom_node_point() + 
    labs(caption = paste0('Layout: ', layout))
  CairoPNG(paste0("ggkar_", layout, ".png"), 800, 800)
  show(gg)
  dev.off()
})
