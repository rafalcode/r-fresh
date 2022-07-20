#!/usr/bin/env Rscript
library(gplots)
library(Cairo)
library(RColorBrewer)

pal <- brewer.pal(11, "RdBu")
pal3 <- pal[c(2,7,11)]
col <- colorRampPalette(pal3)(64)
data(mtcars)

x <- as.matrix(mtcars)
# introduced, early, but, necessary?
# tried, but no reall effect at least in first few examples.
rc <- rainbow(nrow(x), start=0, end=.3)
cc <- rainbow(ncol(x), start=0, end=.3)
##
## demonstrate the effect of row and column dendrogram options
##
CairoPNG("aff.png", 800, 800)
heatmap.2(x, col=col)
dev.off()
# above is not great, data is un normalized?
## default - dendrogram plotted and reordering done.
# heatmap.2(x, dendrogram="none") ## no dendrogram plotted, but reordering done.
# heatmap.2(x, dendrogram="row") ## row dendrogram plotted and row reordering done.
# heatmap.2(x, dendrogram="col") ## col dendrogram plotted and col reordering done.
# heatmap.2(x, keysize=2)
# 
# ## default - dendrogram plotted and reordering done.
# 
# heatmap.2(x, Rowv=FALSE, dendrogram="both") ## generates a warning!
# heatmap.2(x, Rowv=NULL, dendrogram="both") ## generates a warning!
# heatmap.2(x, Colv=FALSE, dendrogram="both") ## generates a warning!
# ## Reorder dendrogram by branch means rather than sums
# heatmap.2(x, reorderfun=function(d, w) reorder(d, w, agglo.FUN = mean) )
# ## plot a sub-cluster using the same color coding as for the full heatmap
# full <- heatmap.2(x)
# heatmap.2(x, Colv=full$colDendrogram[[2]], breaks=full$breaks) # column subset
# heatmap.2(x, Rowv=full$rowDendrogram[[1]], breaks=full$breaks) # row subset
# heatmap.2(x, Colv=full$colDendrogram[[2]],
# Rowv=full$rowDendrogram[[1]], breaks=full$breaks) # both
