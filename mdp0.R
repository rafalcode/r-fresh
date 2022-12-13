#!/usr/bin/env Rscript
# the archetypal example from limma's MDSplot()
# this can serve as simple
library(Cairo)
library(limma)

# Simulate gene expression data for 1000 probes and 6 microarrays.
# Samples are in two groups
# First 50 probes are differentially expressed in second group
rchiv <- rchisq(1000,df=4) # 1000 values
sd <- 0.3*sqrt(4/rchiv) # sure enough, 1000 values too.
x <- matrix(rnorm(1000*6,sd=sd),1000,6)
stop("o")
# so in this way, each gene has a seprate variance, but it's the same for each sample
rownames(x) <- paste("Gene",1:1000)

# first 50 genes (i.e. very few) for the last three samples (Grp2 I suppose) get perturbed.
x[1:50,4:6] <- x[1:50,4:6] + 2 # 2 is just over the max value of sd, ensuring a clear differnce.
# without labels, indexes of samples are plotted.
# mds <- plotMDS(x,  col=c(rep("black",3), rep("red",3)) )
# forget that, go with labels.
# Grp1 in black, grp in red.
CairoPNG("mdp0.png", 800, 800)
# or labels can be provided, here group indicators:
# md <- plotMDS(mds,  col=c(rep("black",3), rep("red",3)), labels= c(rep("Grp1",3), rep("Grp2",3)))
md <- plotMDS(x,  col=c(rep("black",3), rep("red",3)), labels= c(rep("Grp1",3), rep("Grp2",3)))
dev.off()
