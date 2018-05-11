#!/usr/bin/env Rscript
# about named vectors

nv <- 11:15

object.size(nv)
names(nv) <- c("Num1", "Sum1", "Lum1", "Dum1", "Rum1")

# to see named vectors without the names, use unname() ... note the two n's.

object.size(nv)
object.size(names(nv))

# You're able to search on the names, but always using the names() function
grep("[NL].+", names(nv))
nv[sort(names(nv))]
