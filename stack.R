#!/usr/bin/env Rscript
library(stats)

# so this has PlantGrowth dataset available immediately
# it's 30 obs of two variables.
# dead simmple and also the question that comesout is 
# to what extent variable 1 (wieight) is dependent on variable2 (group - also categorical).
# so it's an archetycal anaysis

f <- stats::formula(PlantGrowth)
# this function will take the first column and set it up for dependence-detection on the secnd.
# it looks like string "weight ~ group" but actually R classes it as a "formula" object.

pg <- unstack(PlantGrowth)
# so this separates the obs on the basis of the categorical variable. The are three levels and
# these become columns of a dta aframe.

# stack(pg)  # now put it back together
# stack(pg, select = -ctrl) # omit one vector
