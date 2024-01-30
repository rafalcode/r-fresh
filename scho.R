#!/usr/bin/env Rscript
# this script does what?
# ref. https://mr.schochastics.net/material/netVizR/
library(igraph)
library(ggraph)
library(Cairo)
library(networkdata)

data("got")

# gotu <- upgrade_graph(got)

gotS1 <- got[[1]]
# CairoPNG("fname.png", 800, 800)

# define a custom color palette
# mycols <- c("#1A5878", "#C44237", "#AD8941", "#E99093", "#50594B", "#8968CD", "#9ACD32")
# mycols <- c("Dark Olive Green", "Cornflower Blue", "Aquamarine", "Light Coral", "Navajo White", "Rosy Brown", "Plum")
mycols <- c("Aquamarine", "Cornflower Blue", "Dark Olive Green", "Light Coral", "Navajo White", "Rosy Brown", "Plum")
# I changed the order and there was no problem!
# it did not obey alphabetic
# compute a clustering for node colors
V(gotS1)$clu <- as.character(membership(cluster_louvain(gotS1)))

# compute degree as node size
V(gotS1)$size <- degree(gotS1)

CairoPNG("scho22.png", 800, 800)
gg <- ggraph(gotS1, layout = "stress") +
  geom_edge_link0(aes(edge_linewidth = weight), edge_colour = "grey66") +
  geom_node_point(aes(fill = clu, size = size), shape = 21) +
  geom_node_text(aes(filter = size >= 26, label = name), family = "serif") +
  scale_fill_manual(values = mycols) +
  scale_edge_width(range = c(0.2, 3)) +
  scale_size(range = c(1, 6)) +
  theme_graph() +
  theme(legend.position = "none")
show(gg)
dev.off()
