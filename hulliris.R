#!/usr/bin/env Rscript
# this script does what?
library(ggforce)
library(Cairo)

#setup your ggplot
ggp <- ggplot(iris, aes(Petal.Length, Petal.Width)) +
      geom_mark_hull(aes(fill = Species, label = Species)) +
        geom_point()

#note how the ggplot can be saved as an object 
# and then rendered in a separate paragraph

# render above setup ggp.
CairoPNG("hulliris.png", 800, 800)
show(ggp)
dev.off()
