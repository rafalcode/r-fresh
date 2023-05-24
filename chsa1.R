#!/usr/bin/env Rscript
# the horror 
library(ComplexHeatmap)
library(Cairo)

# trying just an annot.

ha <- HeatmapAnnotation(bar = sample(letters[1:3], 10, replace = TRUE),
    col = list(bar = c("a" = "red", "b" = "green", "c" = "blue")))

colevs = c("a" = "red", "b" = "green", "c" = "blue")

lets <- sample(letters[1:3], 10, replace = TRUE)
letsdf <- data.frame(letters=lets)

ha2 <- HeatmapAnnotation(df=letsdf, col=colevs)

CairoPNG("chsa1.png", 800, 800)
draw(ha2)
dev.off()
