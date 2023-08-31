#!/usr/bin/env Rscript
# this script does what?

mat <- matrix(c(20, 30, 10, 40), byrow=T, nrow=2)
colnames(mat) <- c("Smoking", "Non-Smoking")
rownames(mat) <- c("Male", "Female")

# 
# Male             20              30
# Female           10              40
#

print(fisher.test(mat))
