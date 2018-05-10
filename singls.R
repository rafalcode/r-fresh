#!/usr/bin/env Rscript
# catching singletons
# this is taken from a StOv post, forget which now.
# if anything, can be used to illustrate (albeit not completely) the expand.grid() func.
# Often you have columns that have a unique value which doesn't appear in oter samples
# and you're interested in these unique occurences, or particularly not intersted in them
# let's 
M <- matrix(c(0,1,0,1,0,0,1,0,0,1,0,0,0,1,0,0,0,0,0,0,1,0,0,0,0,0,0,0,1,1,0,0,0,0,1,0,0,0,0,0,1,1), nrow=6)
# So this matrix has rows which register for only one column, and others that register for several.
# so that is how we specifiy singletons in this case, and it can be seen that rows 1 and 3 are the two singletons.
#      [,1] [,2] [,3] [,4] [,5] [,6] [,7]
# [1,]    0    1    0    0    0    0    0
# [2,]    1    0    1    0    0    0    0
# [3,]    0    0    0    1    0    0    0
# [4,]    1    1    0    0    0    0    0
# [5,]    0    0    0    0    1    1    1
# [6,]    0    0    0    0    1    0    1

# find rows which sum to 1
singRows <- which(rowSums(M) == 1)
# find cols which sum to 1
singCols <- which(colSums(M) == 1)

# expand grid is a way unrolling nested loops by preassinging memory for all the combinations
singCombinations <- expand.grid(singRows, singCols)
singCombinations
# singPairRows <- singRows[sapply(singRows, function(singRow) sum(M[,M[singRow,] == 1]) == 1)]
