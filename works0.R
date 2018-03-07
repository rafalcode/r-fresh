#!/usr/bin/env Rscript
# the base library is loaded by default, but I suppose you could be explicit about it:
# library(base)
## be careful with the format: most things in R are floats 
## only integer-valued reals get coerced to integer.
# so this is the way to printf in R.
# needed because of its consistent habit of putting [1] always.
s<-sprintf("%s is %04i feet tall\n", "Sven", 7)
# in R cat will trhow the string as is ... no newlines added
cat(s) 
