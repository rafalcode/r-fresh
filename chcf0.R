#!/usr/bin/env Rscript
# these are scripts tryign to take advantge of the 
# the cell_fun facility in ComplexHEatmap
suppressPackageStartupMessages(library(ComplexHeatmap))
library(Cairo)

# Example data
mat <- matrix(rnorm(100), 10, 10)

# Create heatmap with cell labels
CairoPNG("chcf00.png", 800, 800)
Heatmap(mat,
        name = "Matrix",
        cell_fun = function(j, i, x, y, width, height, fill) {
          grid.text(sprintf("%.2f", mat[i, j]), x, y, gp = gpar(fontsize = 10))
        })
dev.off()

CairoPNG("chcf01.png", 800, 800)
Heatmap(mat,
    cell_fun = function(j, i, x, y, width, height, fill) {
      grid.text(sprintf("%.2f", mat[i, j]), x, y,
                gp = gpar(fontsize = 10,
                          col = ifelse(mat[i, j] > 0, "black", "white"), # changes text colour is value is -ve
                          # fontfamily = "NimbusSans", # NimbusSans no effect
                          fontface = "bold"))
    })
dev.off()
