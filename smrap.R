#!/usr/bin/env Rscript
# this script does what?
# exploring Brad Duthie stuff
# to wit: https://stirlingcodingclub.github.io/simulating_data/index.html
library(ggplot2)
library(Cairo)


# this is your typical correlation sinulation. Firs we start with the predictor variable x1
N   <- 10000
rho <- 0.99
# x1  <- rnorm(n = N, mean = 0, sd = 1)
x1  <- rnorm(N, 0.5, 0.3) # this is the same, as above are defaults.

# OK the depedent variable, this is a stdanrd way to do it.
# x2  <- (rho * x1) + sqrt(1 - rho*rho) * rnorm(n = N, mean = 0, sd = 1);
# amounts to:
rho2 <- rho*rho # come sto .09 in this case
sq1mr <- sqrt(1-rho2) # square root of 1 minus rho squared, comes to .95 in this case
# x2  <- (rho * x1) + sqrt(1 - rho*rho) * rnorm(N)
x2  <- (rho * x1) + sq1mr*rnorm(N, .5, .3)

# what you get here is an x1 of .3 suprise surprise with incredibly low p-value.
# of course you could say, of course, the normal distri "plays into its hands"
# HOWEVER, note how small the rsquared is!
# you can raise it, by bringing rho up to .8, this will give you an rqs of .635 (not that adjusted is the exact same!) 

# Cairo image template
# CairoPNG("fname.png", 800, 800)
# put plot command here
# dev.off()
Corr<-cor(x2, x1, method = "spearman")
CairoPNG("smrapssca.png", 800, 800)
smoothScatter(x1, x2, transformation = function(x) x ^ 0.4,
                colramp = colorRampPalette(c("#000099", "#00FEFF", "#45FE4F",
                                             "#FCFF00", "#FF9400", "#FF3100")),
                xlab = "x1",
                ylab = "x2", xlim=0:1, ylim=0:1,
                font.lab=2, cex.lab=1.4, cex.axis=1.2, cex.main=1.5, main="Paired Sample 496 DNA Methylation Correlation PostTreatment vs. PreTreatment")
  # abline(lm(preresplot1$avg~preresplot$avg-1))
  abline(lm(x2 ~ x1-1))
  abline(a=0, b=1, lty=2)
  text(bquote(r^2 == .(format(Corr, digits=5))), x=0.2,y=0.9, col = "white", cex = 3)
  dev.off()
