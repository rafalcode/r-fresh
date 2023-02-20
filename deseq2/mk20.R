#!/usr/bin/env Rscript
# using DESEq Make example dataset, which SURELY must be a decent simualtor.
# but probably isn't. Howeve ryou would need to refer to it and at least try it.
# however it's docmumented in a lowlevel style
# here I make things more explicit.

# for one, what are betas? Well, they are most probably, the coefficients of the linmod.
# of  course, the interept is betazero, 

# Oh yes, very definitely. The basic b_0 +b_1 * condB
# is very much the departure.
# makes sense too.

library(DESeq2)
library(Cairo)

# here's the canonical call, which al possible option
# n=ngenes (rows)
# m=nsamples (columns)
# betaSD actually by default (and this is risible) ther is no diff exp! Wow it takes some guts to set that as default.

# the first two variables are obvious.
ngenes <- 16
nsamples <- 8
hnsa <- nsamples/2 # half nsamples (for when number of conditions =2)
# betaSD they say is the standard deviation for non-intercept betas, i.e. beta ~ N(0,betaSD)
# he admits intercept is a "beta" betazero.
# and the betas are set via that draw from the normal with that set SD
bSD <- 2 # something like 2 wil give you say 36 genes being over 4 or -4 diffexp.
bSD2 <- 4 # patient
iM <- 4 #intercept meangenes
iSD <- 2
ourfunc <- function(x){4/x+.1} # for dispMeanRel: whatever that means.

# size factors ... ok I normally wouldn't bother with this.
sF <- rep(1, nsamples)
# probably just a way of making the sample look different in the coutn matrix.

# so here's some "leading code"
b0 <- rnorm(ngenes, iM, iSD)
b1 <- rnorm(ngenes, 0, bSD) # this actually the b1 coefficient
b2 <- rnorm(ngenes, 0, bSD2) # this actually the b1 coefficient
be <- cbind(b0, b1)
# so out of this, you can guess we have betazero (the intercept) and beta1, the coeffcient to the chnage of conda to condb.
dispvec <- 2^(be[, 1])
dispr <- ourfunc(dispvec)

# colDats <- DataFrame(cond = factor(rep(c("A", "B"), times = hnsa)), pat=factor(rep(c("P1", "P2"), each=hnsa)))
colDats <- DataFrame(cond = factor(rep(c("A", "B"), times = hnsa)))
rownames(colDats) <- paste0("S", 1:nsamples)
# x <- model.matrix( ~ colDats$cond + colDats$pat)
x <- model.matrix( ~ colDats$cond)

# the following is actually the linear quation
# x= b_0 + b_1 * CondB , so the intercept and second term are "mixed"
xbe <- x %*% t(be) # just formalises the beta (coefficients) matrix 
# Notes: ngenes are columns here ... expression of genes are different, but right now
# it's the saem for all samples. Becareful ... it's easy to forget genes are columns here.

mu <- t(2^xbe * sF) # magnifies the beta as powers of two ... it'll be dragged back with log2 later.
couMat <- matrix(rnbinom(nsamples * ngenes, mu = mu, size = 1/dispr), ncol = nsamples)
# OK what do our rnbinom docs say?
# mu is average number of failues to get size number of successes (which is usuall the set constraint)
# as you can also see, it integerises everything too.
# Unusualy though size is only ngenes long .. 
dds <- DESeqDataSetFromMatrix(countData = couMat, colData=colDats, design = ~ cond)

# To be able to use the PCA plot function. we need to transform via:
# vsd <- varianceStabilizingTransformation(dds)
# no I don't think so, but the rlog() transformation might workq

rld <- rlog(dds)
# CairoPNG("mk2pca.png", 800, 800)
# plotPCA(rld)
# dev.off()
# dds <- DESeq(dds)
# res <- results(dds)
