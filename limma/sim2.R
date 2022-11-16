#!/usr/bin/env Rscript
# simualte's from the ref man of limma
# Simulate a paired experiment with incomplete blocks
# too random, no artifical diffexpr setting.
ngenes <- 50
nsamples <- 12

Block <- c(1,1,2,2,3,3,4,4,5,6,7,8) # note these really are labels: sz 12
Treat <- factor(rep(c(1,2), nsamples/2)) # also labels of course: sz 12
design <- model.matrix(~Treat)
y <- matrix(rnorm(ngenes*nsamples),ngenes,nsamples)
rownames(y) <- paste0("Gene",1:ngenes)

# Estimate the within-block correlation
dupcor <- duplicateCorrelation(y,design,block=Block)
# dupcor$consensus.correlation
# Estimate the treatment effect using both complete and incomplete blocks
fit <- lmFit(y,design,block=Block,correlation=dupcor$consensus)
fit2 <- eBayes(fit)
tt<- topTable(fit2,coef=2)
