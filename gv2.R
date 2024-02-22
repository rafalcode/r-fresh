#!/usr/bin/env Rscript
# this script does what?

library(Gviz)
library(Cairo)

## This is a reference class therefore we show below
## an example from AnnotationTrack
## An empty object
# AnnotationTrack()

## Construct from individual arguments
st <- c(2000000, 2070000, 2100000, 2160000)
ed <- c(2050000, 2130000, 2150000, 2170000)
str <- c("-", "+", "-", "-")
gr <- c("Group1", "Group2", "Group1", "Group3")
annTrack <- AnnotationTrack(start = st, end = ed, strand = str, chromosome = 7,
                            genome = "hg19", feature = "test", group = gr,
                            id = paste("annTrack item", 1:4),
                            name = "generic annotation", stacking = "squish")
## Plotting
CairoPNG("gv20.png", 800, 800)
plotTracks(annTrack)
dev.off()

## Stacking
# it appears here that it can only be done for subtracks within an overall track.
stacking(annTrack)
stacking(annTrack) <- "dense"
CairoPNG("gv21.png", 800, 800)
plotTracks(annTrack)
dev.off()
