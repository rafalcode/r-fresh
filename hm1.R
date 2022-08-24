#!/usr/bin/env Rscript
# exercise with R's basic heatmap function
# datatset used: mtcars.
# it's certainly capable of alot.
library(Cairo)

# for heatmaps, matrices are the best
x <- as.matrix(mtcars)
rc <- rainbow(nrow(x), start = 0, end = .3) # rainbow() from base R
cc <- rainbow(ncol(x), start = 0, end = .3)

CairoPNG("hm10.png", 800, 800)
hv <- heatmap(x, col = cm.colors(256), scale = "column",
              RowSideColors = rc, ColSideColors = cc, margins = c(5,10),
              xlab = "specification variables", ylab = "Car Models",
              main = "heatmap(<Mtcars data>, ..., scale = \"column\")")
dev.off()

utils::str(hv) # the two re-ordering index vectors
## no column dendrogram (nor reordering) at all:
CairoPNG("hm11.png", 800, 800)
heatmap(x, Colv = NA, col = cm.colors(256), scale = "column",
        RowSideColors = rc, margins = c(5,10),
        xlab = "specification variables", ylab = "Car Models",
        main = "heatmap(<Mtcars data>, ..., scale = \"column\")")
dev.off()

## "no nothing"
CairoPNG("hm12.png", 800, 800)
heatmap(x, Rowv = NA, Colv = NA, scale = "column",
        main = "heatmap(*, NA, NA) ~= image(t(x))")
dev.off()

# round(Ca <- cor(attitude), 2)
# round(Ca <- cor(attitude), 2)
# round(Ca <- cor(attitude), 2)
# symnum(Ca) # simple graphic
# heatmap(Ca, symm = TRUE, margins = c(6,6)) # with reorder()
# heatmap(Ca, Rowv = FALSE, symm = TRUE, margins = c(6,6)) # _NO_ reorder()
# 
# ## slightly artificial with color bar, without and with ordering:
# cc <- rainbow(nrow(Ca))
# heatmap(Ca, Rowv = FALSE, symm = TRUE, RowSideColors = cc, ColSideColors = cc,
# margins = c(6,6))
# heatmap(Ca,symm = TRUE, RowSideColors = cc, ColSideColors = cc,
# margins = c(6,6))
# ## For variable clustering, rather use distance based on cor():
# symnum( cU <- cor(USJudgeRatings) )
# hU <- heatmap(cU, Rowv = FALSE, symm = TRUE, col = topo.colors(16),
#              distfun = function(c) as.dist(1 - c), keep.dendro = TRUE)
# ## The Correlation matrix with same reordering:
# round(100 * cU[hU[[1]], hU[[2]]])
# ## The column dendrogram:
# utils::str(hU$Colv)
# # Cairo image template
# # put plot command here
# dev.off()
