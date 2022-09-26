#!/usr/bin/env Rscript
# a heatmap using pheatmap which is faous fo rits annotation bars.
# it's Kevin Blighe
# ref. https://www.biostars.org/p/351551/#351583
library(pheatmap)
library(Cairo)

# Create random data
data <- replicate(20, rnorm(50))
rownames(data) <- paste("Gene", c(1:nrow(data)))
colnames(data) <- paste("Sample", c(1:ncol(data)))

# Create col and row metadata
metadata <- data.frame(c(rep("case", ncol(data)/2), rep("control", ncol(data)/2)),
      c(rep("cond1", ncol(data)/4), rep("cond2", ncol(data)/4), rep("cond3", ncol(data)/4), rep("cond4", ncol(data)/4)),
      row.names=colnames(data))
colnames(metadata) <- c("casecontrol","condition")

metadata_gene <- data.frame(c(rep("Tcell", nrow(data)/2), rep("Bcell", nrow(data)/2)),
      row.names=rownames(data))
colnames(metadata_gene) <- c("Cell")

# create the heatmap
CairoPNG("ph0.png", 800, 800)
pheatmap(data, show_rownames=F, cluster_cols=T, cluster_rows=T, scale="row",
      cex=1, clustering_distance_rows="euclidean", cex=1,
      clustering_distance_cols="euclidean", clustering_method="complete", border_color=FALSE,
      annotation_col=metadata, annotation_row=metadata_gene)
dev.off()
