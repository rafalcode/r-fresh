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

# Clustering
# mat <- t(mat) # To flip the figure, if needed
rowDendrogram <- as.dendrogram(hclust(dist(mat)))
rowOrder <- order.dendrogram(rowDendrogram)
columnDendrogram <- as.dendrogram(hclust(dist(t(mat))))
columnOrder <- order.dendrogram(columnDendrogram)

# Plot
lattice.options(axis.padding=list(factor = 0.5)) # Sets padding between heatmap and box

# levelplot is one lattice's functions, not sure what's special abou tit.
p1 <- levelplot(mat[rowOrder, columnOrder],
                aspect = "fill", 
                scales = list(x = list(rot = 90), tck = c(1,0)),
                # Specify key color ramp
                colorkey = list(space = "left", col = colorRampPalette(c("blue", "red", "yellow"), space = "Lab")(100), height = 0.75),
                # Specify plot color ramp (same as key)
                col.regions = colorRampPalette(c("blue", "red", "yellow"), space = "Lab")(100),
                xlab = NULL, ylab = NULL, pretty = TRUE, 
                main = "Byproduct ratios in treated water samples",
                # Set range and number of breakpoints in intensity data,
                # length.out should match the number of number of colors
                # specified for the ke and plot.
                at = seq(0, 4500, length.out = 100),
                # Set dendogram position
                par.settings = list(layout.widths = list(axis.key.padding = -1.5), layout.heights = list(key.axis.padding = -1.5)),
                # Draw dendograms
                legend = list(right = list(fun = dendrogramGrob, args = list(x = columnDendrogram, side = "right", size = 5)),                      
                       top = list(fun = dendrogramGrob, args = list(x = rowDendrogram, side = "top", size = 5)))
)
CairoPNG("hml0.png", 800, 800)
show(p1)
dev.off()

# with lattice I think you're able to work on the components separately
