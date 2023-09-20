#!/usr/bin/env Rscript
# this script does what?
# fgsea .. time series RNA
library(Cairo)
library(fgsea)
library(data.table)
library(ggplot2)

data(examplePathways)
# you get:  a list of 1457 lists: it's all the reactome for mouse I think.

data(exampleRanks)
# you get: actually a simple floatin gpoint vector 12000 long.
# names look entrez style. Values -63.33 to +53.28 starting ta lowest.

set.seed(42)

fgseaRes <- fgsea(pathways = examplePathways, 
                  stats = exampleRanks,
                  minSize=15, maxSize=500)
# you get:
# a struct which has an $NES
# and a list of 586 lists calle leadingEdge

