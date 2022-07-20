#!/usr/bin/env Rscript
library(Cairo)
# heatmaps with lattice
library(lattice)
library(latticeExtra)

## ---- fig.height=6.5, fig.width=6----------------------------------------
# Read data old style, I comment out and show the better row.names=1 method.
# x <- read.csv("hml_exampda1.csv", stringsAsFactors = FALSE)
# Set up matrix by moving sample codes from first column to row names
# sampleCodes <- x[, 1] # Grab sample codes
# mat <- data.matrix(x[, -1]) # Remove first column
# rownames(mat) <- sampleCodes # Assign as row names
mat <- as.matrix(read.csv("hml_exampda1.csv", row.names=1)) # Rv4 way: above 4 lines into 1.

# Clustering and preparation
# mat <- t(mat) # To flip the figure, if needed
rowDendrogram <- as.dendrogram(hclust(dist(mat)))
rowOrder <- order.dendrogram(rowDendrogram)
columnDendrogram <- as.dendrogram(hclust(dist(t(mat))))
columnOrder <- order.dendrogram(columnDendrogram)

# Plot
lattice.options(axis.padding=list(factor = 0.5)) # Sets padding between heatmap and box

p2 <- levelplot(mat[rowOrder, columnOrder],
          aspect = "fill", scales = list(x = list(rot = 90), tck = c(1,0)),
          colorkey = list(space = "left", col = colorRampPalette(c("blue", "red", "yellow"), space = "Lab")(10), height = 0.75),
          col.regions = colorRampPalette(c("blue", "red", "yellow"), space = "Lab")(10),
          # diff with hml0.R, just (10) here insteadof (100)
          # next, cuts also set to 10.
          xlab = NULL, ylab = NULL, cuts = 10, pretty = TRUE,
          main = "Byproduct ratios in treated water samples",
          par.settings = list(layout.widths = list(axis.key.padding = -1.5), layout.heights = list(key.axis.padding = -1.5)),
          legend = list(right = list(fun = dendrogramGrob, args = list(x = columnDendrogram, side = "right", size = 5)),                      
                 top = list(fun = dendrogramGrob, args = list(x = rowDendrogram, side = "top", size = 5)))
)
CairoPNG("hml1.png", 800, 800)
show(p2)
dev.off()
