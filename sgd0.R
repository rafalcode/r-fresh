#!/usr/bin/env Rscript
# G for Gerard, G for garbled.
library(seqgendiff)
library(Cairo)

## Simulate data from given matrix of counts
## In practice, you would obtain Y from a real dataset, not simulate it.
set.seed(1)
nsamp <- 10
ngene <- 1000
Y <- matrix(stats::rpois(nsamp * ngene, lambda = 50), nrow = ngene)


thinout <- thin_2group(mat=Y, prop_null=0.9, signal_fun=stats::rexp, signal_params=list(rate = 0.5))
# prop_null clearly is proportional of nulls or zero.
# very cavalier with the design matrix too, it's that signal params,
# basically 50% in COndA, 50% in CindB distributed randomly! Talk about not giving a shit!

# what does null mean? Annoying terminology by Mr. Gerard
# it's clear it's zero
# and of coure, itcan't mean any other thing except zero count genes!
# Check 90 percent of genes are null:
# mean(abs(thinout$coef) < 10^-6)
# or thinout$coefmat.

## Check the estimates of the log2-fold change
Ynew <- log2(t(thinout$mat + 0.5)) # 0.5 to protect against log2(0) of course.

X    <- thinout$designmat
Bhat <- coef(lm(Ynew ~ X))["X", ]

CairoPNG("sgd0.png", 800, 800)
plot(thinout$coefmat, Bhat, xlab = "Actual log2FC", ylab = "Estimated log2FC")
abline(0, 1, col = 2, lwd = 2)
dev.off()
