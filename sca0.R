#!/usr/bin/env Rscript
# so what exactly does scale() do, and how can it do it row wise
# considering its help file talks about doing it 
library(Cairo)

set.seed(100)
a<- sort(rnorm(8,1,0.5))
b<- sort(rnorm(8,2,1))
c<- sort(rnorm(8,4,2))
d<- sort(rnorm(8,9,3))
ma <- matrix(c(a,b,c,d), ncol=4)
# Cairo image template
# CairoPNG("fname.png", 800, 800)
# put plot command here
# dev.off()
