#!/usr/bin/env Rscript
# this script does what?
library(ggplot2)
library(Cairo)
library(grid)

CairoPNG("gvp.png", 400, 400)

stripVP <- viewport(y=1, height=unit(1, "lines"), just="top", name="stripvp")
panelVP <- viewport(y=0, height=unit(1, "npc") - unit(1, "lines"), just="bottom", name="panelvp")
pushViewport(stripVP)
grid.rect(gp=gpar(fill="grey80"), name="striprect")
upViewport()
pushViewport(panelVP)
grid.rect(name="panelrect")
upViewport()

downViewport("stripvp")
grid.text("strip text", name="striptext")
upViewport()

grid.edit("striptext", label="modified text", gp=gpar(col="darkgreen", fontface="italic", fontfamily="serif"))
grid.ls(grid.force()) # "col_annotation" looks like it's the one to edit
dev.off()
