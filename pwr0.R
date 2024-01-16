#!/usr/bin/env Rscript
# the example from the pwr vignette

library(pwr)

efsz <- ES.h(p1 = 0.75, p2 = 0.50) # p is for proportion 

res <- pwr.p.test(h=efsz, sig.level=.05, power=.8, alternative="greater")

# efsz, the effectsize is cirtical, as themanual says:
# set p1 to .65 ... i.e. make effect size smaller, What do you get for n?
# About 85 coin flips. Detecting smaller effects require larger sample sizes!

# that's the spanner in the Power analysis works!
