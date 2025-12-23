# using complexHeatmap for triangular heatmaps.

library(ComplexHeatmap)
library(Cairo)

m <- cor(matrix(rnorm(100), 10))
CairoPNG("cohm20.png", 800, 800)
Heatmap(m)
dev.off()

# By observing the heatmap, the simplest way is to compare the row index and column index of the heatmap. Note here since the heatmap rows and columns are already reordered by clusering, we need to reorder the matrix before sending to heatmap, and in the heatmap, no reordering should be applied.

od =  hclust(dist(m))$order
m2 = m[od, od]

CairoPNG("cohm21.png", 800, 800)
Heatmap(m2, rect_gp = gpar(type = "none"), 
    cluster_rows = FALSE, cluster_columns = FALSE,
    # now the following realy is tricky to understand, some variables not even declared.
    cell_fun = function(j, i, x, y, w, h, fill) {
        if(i >= j) {
            grid.rect(x, y, w, h, gp = gpar(fill = fill, col = fill))
        }
    })
dev.off()
