#!/usr/bin/env Rscript
library(beeswarm)

myz = rnorm(100000,0,1) 

mydata = sample(myz, 20) 
mygroup = c(rep('A', 10), rep('B', 10))

mydata[11:20] = mydata[11:20] + 1

beeswarm(mydata ~ mygroup)
