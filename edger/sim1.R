#!/usr/bin/env Rscript
# n the refman for edgeR these are only two mentions of simulation and this is one of them, for the MDS plot.
# which is a multi-dimensional scaling plot.
library(edgeR)
library(Cairo)

# calctpm1s: "calculate the tpm for 1 sample" (.e. 1 column vector...  we'll need apply() func for a matrix)
calctpm1s <- function(counts, lengths)
{
  rate <- counts / lengths
  rate / sum(rate) * 1e6
}

# Simulate DGE data for 1000 genes and 6 samples.
# Samples are in two groups
# First 200 genes are differentially expressed in second group
ngenes <- 1000
lgenes <- round(1000*runif(ngenes)+1) # length of each of the genes
nlib <- 6 # aka. nsamps
counts <- matrix(rnbinom(ngenes*nlib, size=1/10, mu=20),ngenes,nlib)
rownames(counts) <- paste("gene",1:ngenes, sep=".")

# two groups which behave differently
group <- gl(2,3,labels=c("Grp1","Grp2")) # gl() generate factor levels
# In group 2, first 200 genes will be differentially expressed purposely.
counts[1:200,group=="Grp2"] <- counts[1:200,group=="Grp2"] + 10
# lets check out tpms at this stage:
tpms <- apply(counts, 2, function(x) calctpm1s(x, lgenes))

y <- DGEList(counts,group=group)
y2 <- calcNormFactors(y) # of course this does nto take into account gene lengths.
# without labels, indexes of samples are plotted.
col <- as.numeric(group)

CairoPNG("sim1a.png", 800, 800)
mds <- plotMDS(y2, top=200, col=col)
# show(mds) # will output stuff to console.
dev.off()

# or labels can be provided, here group indicators:
CairoPNG("sim1b.png", 800, 800)
plotMDS(mds, col=col, labels=group)
dev.off()
