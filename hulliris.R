#!/usr/bin/env Rscript
# this is from ggforce's
# the "force" name ain't great .. it's a sort of accelerator
# of ggplot whith a diverse set of functions

# the hull ringing feautre is good... -ish, but
it can;t ope with fairly difficult ones at all.

# ANd in fact this is very apparent in the example of iris
# a bunch of the versicolors and stuck right in virginica territory
# and the hull can't go in and out grabbing each one, so this is when ovberalps occur.

# I don;t like the hull labels and they are a devil to get rid of
# actually the trick is to get them (i.e. option label=) out of the main ggplot aes() section
library(ggforce)
library(Cairo)

ir <- iris[c(25:30, 85:90, 112:130),]

ggp <- ggplot(ir, aes(Petal.Length, Petal.Width)) +
      # geom_mark_hull(aes(fill = Species, filter = Species != 'versicolor')) +
      # geom_mark_hull(aes(fill = Species)) +
      geom_mark_hull(aes(fill = Species)) +
        geom_point()
#setup your ggplot
ggp2 <- ggplot(ir, aes(Petal.Length, Petal.Width)) + 
      # geom_mark_hull(aes(group = Species, label = Species)) +
      geom_mark_hull(aes(group = Species)) +
        geom_point() +
        # geom_text(label=rownames(iris))
        geom_text(label=rownames(ir), hjust=1.5, vjust=1.5)

#note how the ggplot can be saved as an object 
# and then rendered in a separate paragraph

# render above setup ggp.
CairoPNG("hulliris.png", 800, 800)
show(ggp2)
dev.off()
