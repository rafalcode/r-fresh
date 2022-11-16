#!/usr/bin/env Rscript
# another of the simualtions in the limma ref.
# gives a heatmap ... a bit rough mind you.
library(ggplot2)
library(Cairo)

# Simulate gene expression data for 50 genes and 6 microarrays.
# Samples are in two groups
# First 50 probes are differentially expressed in second group
ngenes <- 50
sd <- 0.3*sqrt(4/rchisq(ngenes,df=4)) #varies between 0.1 and 0.7 mostly in 0.1 to 0.33
x <- matrix(rnorm(ngenes*6,sd=sd),ngenes,6)
rownames(x) <- paste("Gene",1:ngenes)
x <- x + seq(from=0, to=16, length=ngenes)
x[,4:6] <- x[,4:6] + 2 # easy but enitrely unrealistic.

# Cairo image template
CairoPNG("sim1.png", 800, 800)
coolmap(x)
dev.off()
