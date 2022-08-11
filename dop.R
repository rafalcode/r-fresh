#!/usr/bin/env Rscript
# example fo do[parallel in the doparallel_vignette
library(doParallel)

# registering paralell environment ... I think snow works like this 
cl <- makeCluster(4)
registerDoParallel(cl)

x <- iris[which(iris[,5] != "setosa"), c(1,5)]
trials <- 10000
ptime <- system.time({r <- foreach(icount(trials), .combine=cbind) %dopar% {
        ind <- sample(100, 100, replace=TRUE)
        result1 <- glm(x[ind,2]~x[ind,1], family=binomial(logit))
        coefficients(result1)
        }})[3]

# on p79t I got 10.2 secs for 2 core clust
# 5.9 with 4 core cluster.
