# example from the ENmix package

library(ENmix)
library(minfiData)

path <- file.path(find.package("minfiData"),"extdata")
rgSet <- readidat(path = path,recursive = TRUE)
mdat <- getmeth(rgSet)

beta <- getB(mdat,"Illumina")
# with that you get a 485577x6 numeric matrix, so it's the 450k chip with 6 samples.

m <- B2M(beta)

# and back again.
b <- M2B(m)
# yes looks the exact same.
# I know it's an easy mathematical formula. no time to look at that right now.
