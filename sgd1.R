#!/usr/bin/env Rscript
# G for Gerard, G for garbled.
# (later) ha!
library(seqgendiff)
library(Cairo)

# For the special case when your design matrix is just a group indicator (that is, you have two groups of individuals), you can use the function thin_2group(). Let’s generate data from the two-group model where 90% of genes are null and the non-null effects are gamma-distributed.

# rgamma used this time!
thout <- thin_2group(mat = submat, 
                     prop_null     = 0.9, 
                     signal_fun    = stats::rgamma,
                     signal_params = list(shape = 1, rate = 1))

# We can again verify that we thinned appropriately using the voom-limma pipeline:

new_design <- cbind(thout$design_obs, thout$designmat)
# new_design
#>      (Intercept) P1
#> [1,]           1  0
#> [2,]           1  0
#> [3,]           1  0
#> [4,]           1  1
#> [5,]           1  1
#> [6,]           1  1

vout <- limma::voom(counts = thout$mat, design = new_design)
lout <- limma::lmFit(vout)
coefhat <- coef(lout)[, 2, drop = FALSE]
And we can plot the results


CairoPNG(sgd1.png, 800, 800)
oldpar <- par(mar = c(2.5, 2.5, 1, 0) + 0.1, mgp = c(1.5, 0.5, 0))
plot(x    = thout$coefmat, 
     y    = coefhat, 
     xlab = "True Coefficient", 
     ylab = "Estimated Coefficient",
     main = "First Variable",
     pch  = 16)
abline(a   = 0,
       b   = 1,
       lty = 2, 
       col = 2,
       lwd = 2)
dev.off()
