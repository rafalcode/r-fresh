#!/usr/bin/env Rscript
# example from refman ref. generateSyntheticData() function.
# proves you don't need a
library(ggplot2)
library(Cairo)
library(compcodeR)
# library(ape)

## RNA-Seq data
ccda <- generateSyntheticData(dataset = "mydata", n.vars = 100, samples.per.cond = 3, n.diffexp = 50)
# in this case first 3 samples will be first condition
# first 50 of genes will be diffexp.
# is dataset option actually needed?
# ccda <- generateSyntheticData(n.vars = 100, samples.per.cond = 3, n.diffexp = 50)
# answer to that is yes ... you will get an error with out dataset option defined.

# ccda@count.matrix will get you the count matrix, genes (rows) called g1, g2, etc. cols sample1, etc.
## Inter-species RNA-Seq data
# tree <- read.tree(text = "(((A1:0,A2:0,A3:0):1,B1:1):1,((C1:0,C2:0):1.5,(D1:0,D2:0):1.5):0.5);")
# id.species <- factor(c("A", "A", "A", "B", "C", "C", "D", "D"))
# names(id.species) <- tree$tip.label
# mydata.obj <- generateSyntheticData(dataset = "mydata", n.vars = 1000,
#                                      samples.per.cond = 4, n.diffexp = 100,
#                                      tree = tree,
#                                      id.species = id.species,
#                                      lengths.relmeans = "auto",
#                                      lengths.dispersions = "auto")
# Cairo image template
# CairoPNG("fname.png", 800, 800)
# put plot command here
# dev.off()
