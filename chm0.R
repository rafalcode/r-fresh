#!/usr/bin/env Rscript
# this script does what? complexHeatmap start of tute
# although what I'm trying to get it the blue to yellow.
library(ggplot2)
library(Cairo)
library(ComplexHeatmap)
library(circlize)

set.seed(123)
nr1 = 4
nr2 = 8
nr3 = 6
nr = nr1 + nr2 + nr3

nc1 = 6
nc2 = 8
nc3 = 10
nc = nc1 + nc2 + nc3

# kinf oquite a complicated way:
mat = cbind(rbind(matrix(rnorm(nr1*nc1, mean = 1,   sd = 0.5), nr = nr1),
          matrix(rnorm(nr2*nc1, mean = 0,   sd = 0.5), nr = nr2),
          matrix(rnorm(nr3*nc1, mean = 0,   sd = 0.5), nr = nr3)),
    rbind(matrix(rnorm(nr1*nc2, mean = 0,   sd = 0.5), nr = nr1),
          matrix(rnorm(nr2*nc2, mean = 1,   sd = 0.5), nr = nr2),
          matrix(rnorm(nr3*nc2, mean = 0,   sd = 0.5), nr = nr3)),
    rbind(matrix(rnorm(nr1*nc3, mean = 0.5, sd = 0.5), nr = nr1),
          matrix(rnorm(nr2*nc3, mean = 0.5, sd = 0.5), nr = nr2),
          matrix(rnorm(nr3*nc3, mean = 1,   sd = 0.5), nr = nr3))
   )
mat = mat[sample(nr, nr), sample(nc, nc)] # random shuffle rows and columns
rownames(mat) = paste0("row", seq_len(nr))
colnames(mat) = paste0("column", seq_len(nc))

col_fun = colorRamp2(c(-1, 0, 1), c("#0101ed", "#818175", "#fbfb02"))
# yes this seem to be function, you would call it manually this way.
# col_fun(seq(-3, 3))

hm <- ComplexHeatmap::Heatmap(mat, name = "mat", col = col_fun)

# Cairo image template
CairoPNG("chm0.png", 800, 800)
show(hm)
dev.off()
