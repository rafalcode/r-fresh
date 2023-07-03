#!/usr/bin/env Rscript
# you don't acutally need tibbles to use gather()
library(ggplot2)
library(tidyr)
library(Cairo)

df2 <- data.frame(player=c('A', 'B', 'C', 'D'),
                  year1=c(12, 15, 19, 19),
                  year2=c(22, 29, 18, 12),
                  year3=c(17, 17, 22, 25))
# ga <- df2 %>% gather(key="year", value="points", 2:4)
# that works
ga <- df2 %>% gather(year, points, 2:4)
# the same
ga2 <- df2 %>% gather(year, points, year1:year2)
# we see here how year3 is kept as is. That it's simply repeated.

# And note what a misnomer gather is. There is no collecting into 1 entity unless you called the 
# column year and entity. There is no adding of numbers either.
# It's just another hacky way of creating long form

# Of course all this could make sense in the context of the ggplot
# which requires long form

# Cairo image template
# CairoPNG("fname.png", 800, 800)
# put plot command here
# dev.off()
