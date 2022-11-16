#!/usr/bin/env Rscript
# this script does what?
library(limma)
library(Cairo)

#  See also lmFit examples
# Simulate gene expression data,
# 6 microarrays and 100 genes with one gene differentially expressed
set.seed(2016)

sigma2 <- 0.05 / rchisq(100, df=10) * 10
y <- matrix(rnorm(100*6,sd=sqrt(sigma2)),100,6)
design <- cbind(Intercept=1,Group=c(0,0,0,1,1,1))

# perturb first gene on last 3 samples
y[1,4:6] <- y[1,4:6] + 1

fit <- lmFit(y,design)
# Moderated t-statistic
fit <- eBayes(fit)
tt0 <- topTable(fit,coef=2)

# Ordinary t-statistic: but this isn't used?
ordinary.t <- fit$coef[,2] / fit$stdev.unscaled[,2] / fit$sigma

# Treat relative to a 10% fold-change
# use of treat here hard to understand.
tfit <- treat(fit, fc=1.1)
tt2 <- topTreat(tfit,coef=2)
