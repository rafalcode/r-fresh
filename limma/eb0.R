#!/usr/bin/env Rscript
# the ebayes examplesin the limma refman.
library(limma)
library(Cairo)

#  See also lmFit examples
# Simulate gene expression data,
# 6 microarrays and 100 genes with one gene differentially expressed
set.seed(2016)
sigma2 <- 0.05 / rchisq(100, df=10) * 10
y <- matrix(rnorm(100*6,sd=sqrt(sigma2)),100,6)
design <- cbind(Intercept=1,Group=c(0,0,0,1,1,1))
y[1,4:6] <- y[1,4:6] + 1
fit <- lmFit(y,design)

# Moderated t-statistic
fit2 <- eBayes(fit)
topTable(fit2,coef=2)

# Question: how does fit differ from fit2?
# Well both are MArrayLM objects:
# but can call each of their components by names if you use double squarebrackets or $ referencers.
# names() will get you those names, fit just has 12: 
# coefficients rank assign qr df.residual sigma cov.coefficients stdev.unscaled pivot Amean method design
# fit2 has an additional 11:
# df.prior s2.prior var.prior proportion s2.post t
# df.total p.value lods F F.p.value

# Ordinary t-statistic
ordinary.t <- fit2$coef[,2] / fit2$stdev.unscaled[,2] / fit2$sigma

# Treat relative to a 10% fold-change
tfit <- treat(fit2, fc=1.1)
# tfit very fit2
topTreat(tfit,coef=2)
