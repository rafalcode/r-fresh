#!/usr/bin/env Rscript
# I'm trying to understand this grob thing and how to apply it
library(Cairo)
library(grid)

# CairoPNG("grob0.png", 800, 800)
# grid.rect(gp=gpar(col="red"))
# same thing, but more verbose
grid.draw(applyEdit(rectGrob(), gEdit(gp=gpar(col="red"))))
# Cairo image template
# put plot command here
# dev.off()

# Learnings:
# Yes, 
# grid.rect(gp=gpar(col="red"))
# and
# grid.draw(applyEdit(rectGrob(), gEdit(gp=gpar(col="red"))))
# do prodcue exactly the same thing.
# it does produce and image of white fill
# about 600 or 800 square (a default at work?)
# (ehem .. it's 800. That's what you set your CairoPNG to)
# you can faintly detect red on the outside

