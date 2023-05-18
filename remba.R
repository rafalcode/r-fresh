#!/usr/bin/env Rscript
# the example for removeBatchEffect()
# it's very much throwaway, just use rnorm to simulate.
library(limma)
library(Cairo)

nsamps <- 9
nfeatures <- 10
y <- matrix(rnorm(nfeatures*nsamps),nfeatures,nsamps) # nfeatures will be the number of rows.
batch <- c("A","A","A","B","B","B","C","C","C")
y[,1:3] <- y[,1:3] + 5
# because of this arrangemet remBE will wipe out the diferrent I expect.
# that said there is no bias between B's and C's ... if I put one in, would it affect?

y2 <- removeBatchEffect(y, batch) # CRITICAL, the batch should match the columns!

CairoPNG("remba.png", 800, 800)
par(mfrow=c(1,2))
boxplot(as.data.frame(y),main="Original")
boxplot(as.data.frame(y2),main="Batch corrected")
dev.off()
