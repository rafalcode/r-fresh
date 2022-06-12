#!/usr/bin/env Rscript
# this script does what?
library(ggplot2)
library(Cairo)
library(sva)
library(SmartSVA)

## Methylation M values (CpG by Sample)
Y <- matrix(rnorm(20*1000), 1000, 20)
df <- data.frame(pred=gl(2, 10))
## Determine the number of SVs
Y.r <- t(resid(lm(t(Y) ~ pred, data=df)))
## Add one extra dimension to compensate potential loss of 1 degree of freedom
## in confounded scenarios (very important)
n.sv <- EstDimRMT(Y.r, FALSE)$dim + 1
mod <- model.matrix( ~ pred, df)
svo <- smartsva.cpp(Y, mod, mod0=NULL, n.sv=n.sv)
sv <- sva(Y, mod, n.sv=n.sv)

## Speed comparison to traditional SVA
## Not run:
## Methylation M values (CpG by Sample, 27K by 1,000)
Y <- matrix(rnorm(1000*27000), 27000, 1000)
df <- data.frame(pred=gl(2, 500))
## Determine the number of SVs
Y.r <- t(resid(lm(t(Y) ~ pred, data=df)))
n.sv <- 50
mod <- model.matrix( ~ pred, df)
system.time(sv.obj1 <- smartsva.cpp(Y, mod, mod0=NULL, B=5, alpha = 1, VERBOSE=TRUE, n.sv=n.sv))
system.time(sv.obj2 <- sva(Y, mod, mod0=NULL, B=5, n.sv=n.sv))
## Check if the solutions are the same
# head(sv.obj1$sv)
# head(sv.obj2$s
# actually prtty close here, but not above!
# Cairo image template
# CairoPNG("fname.png", 800, 800)
# put plot command here
# dev.off()
