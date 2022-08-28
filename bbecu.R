#!/usr/bin/env Rscript
# this hard coded for Davd Robo's curves (cu)
# becu is for beta curves.
# it's bbecu.R, the baseball becu.
library(Cairo)

sh1 <- 71 # average made bats
sh2 <- 291 # average at-bats

# lucky starter
made1 <- 4 # bats made
att1 <- 10 # bats attempted.

# true veteran
made2 <- 300 # bats made
att2 <- 1000 # bats attempted.

sq <- seq(0, 1, by=0.005)
bvals <- dbeta(sq, sh1, sh2) # uniform beta vals (the x is uniform, not the y!)

bv1 <- dbeta(sq, sh1+made1, sh2+att1-made1)
bv2 <- dbeta(sq, sh2+made2, sh2+att2-made2)

# Beware when looking at this values, these are very thin densities, only a certain narrow
# setion (around 0.27 or so) will have decent values, all other with be miniscule and
# so will look odd.
mm <- cbind(bvals, bv1, bv2)

CairoPNG("bb0.png", 800, 800)
plot(sq, bvals, pch=20, cex=0.5, col="steelblue", type="l")
dev.off()

CairoPNG("bb1.png", 800, 800)
plot(sq, bv1, pch=20, cex=0.5, col="goldenrod", type="l")
dev.off()

CairoPNG("bb2.png", 800, 800)
plot(sq, bv2, pch=20, cex=0.5, col="firebrick", type="l")
dev.off()

CairoPNG("bb3.png", 800, 800)
matplot(mm, type="l")
dev.off()
