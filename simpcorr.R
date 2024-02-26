#!/usr/bin/env Rscript
# this script does what?
# exploring Brad Duthie stuff
# to wit: https://stirlingcodingclub.github.io/simulating_data/index.html
# trying to do smmothscatter plots and work them out.

# I'm sort of getting expected resutls here ...
# though I'm not sure
library(ggplot2)
library(Cairo)


# this is your typical correlation sinulation. Firs we start with the predictor variable x1
rho <- 0.83
x1  <- rnorm(500, mean = 0, sd = 1)
# x1  <- 1:10000

# OK the depedent variable, this is a stdanrd way to do it.
# x2  <- (rho * x1) + sqrt(1 - rho*rho) * rnorm(n = N, mean = 0, sd = 1);
# amounts to:
rho2 <- rho*rho # come sto .09 in this case
sq1mr <- sqrt(1-rho2) # square root of 1 minus rho squared, comes to .95 in this case
# x2  <- (rho * x1) + sqrt(1 - rho*rho) * rnorm(N)
x2  <- (rho * x1) + sq1mr*rnorm(length(x1), .5, .3)

# despite this cofident pre-setting of correlation, it's only being pertrubed by normally distributed errors.


Corr<-cor(x2, x1, method = "spearman")
CairoPNG("smpcorr0.png", 800, 800)
smoothScatter(x1, x2, transformation = function(x) x ^ 0.4,
                colramp = colorRampPalette(c("#000099", "#00FEFF", "#45FE4F",
                                             "#FCFF00", "#FF9400", "#FF3100")),
              bandwidth=5,
                xlab = "x1",
                ylab = "x2", xlim=0:1, ylim=0:1,
                font.lab=2, cex.lab=1.4, cex.axis=1.2, cex.main=1.5, main="Sim. of Correlation")
  # abline(lm(preresplot1$avg~preresplot$avg-1))
  # abline(lm(x2 ~ x1-1))
  # abline(a=0, b=1, lty=2)
  abline(a=0, b=1)
  text(bquote(r^2 == .(format(Corr, digits=5))), x=0.2,y=0.9, col = "white", cex = 3)
  dev.off()
df <- data.frame(x1=x1, x2=x2)
CairoPNG("smpcorr1.png", 800, 800)
gp <- ggplot(df, aes(x=x1, y=x2)) +
    geom_point(alpha=.2, size=3)
show(gp)
dev.off()

# scatter plots colorimetrics: I already hit up on it, they end up being symmetrical
