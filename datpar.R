#!/usr/bin/env Rscript
# this script does what?
# experieting with carte's creat data partition
# get nothing out of this.
library(ggplot2)
library(Cairo)
library(caret)

data(oil) # get fattyAcids and oilType from this (only!)
oilp <- createDataPartition(oilType, 2)

x <- rgamma(50, 3, .5) # 50 values
inA <- createDataPartition(x, list = FALSE) # indices for 26 rows

CairoPNG("datpar.png", 800, 800)
plot(density(x[inA]))
dev.off()
stop("o")
rug(x[inA])

points(density(x[-inA]), type = "l", col = 4)
rug(x[-inA], col = 4)

createResample(oilType, 2)

createFolds(oilType, 10)
createFolds(oilType, 5, FALSE)

createFolds(rnorm(21))

createTimeSlices(1:9, 5, 1, fixedWindow = FALSE)
createTimeSlices(1:9, 5, 1, fixedWindow = TRUE)
createTimeSlices(1:9, 5, 3, fixedWindow = TRUE)
createTimeSlices(1:9, 5, 3, fixedWindow = FALSE)

createTimeSlices(1:15, 5, 3)
createTimeSlices(1:15, 5, 3, skip = 2)
createTimeSlices(1:15, 5, 3, skip = 3)

set.seed(131)
groups <- sort(sample(letters[1:4], size = 20, replace = TRUE))
# table(groups)
folds <- groupKFold(groups)
# lapply(folds, function(x, y) table(y[x]), y = groups
