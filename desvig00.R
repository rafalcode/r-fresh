#!/usr/bin/env Rscript
# http://bioconductor.org/packages/devel/bioc/vignettes/DESeq2/inst/doc/DESeq2.html
# This is the original code, in the vignette Rmd file in the source code.
# desvig0, one zero, is the modified version.
# but echo=F for the RMD so this code doesn't come out.

library(ggplot2)
library(Cairo)
library(gridExtra)

plotit0 <- function(d, title) {
  # ggplot(d, aes(x=cond, y=log2c, group=1)) + geom_point() +
  ggplot(d, aes(x=cond, y=log2c)) + geom_point() +
    stat_summary(fun=mean, geom="line", colour="red", size=0.8) + 
    xlab("condition") + ylab("log2(counts+1)") + ggtitle(title)
}

plotit <- function(d, title) {
  ggplot(d, aes(x=cond, y=log2c, group=geno)) + 
    geom_jitter(size=1.5, position = position_jitter(width=.15)) +
    facet_wrap(~ geno) + 
    stat_summary(fun=mean, geom="line", colour="red", size=0.8) + 
    xlab("condition") + ylab("log2(counts+1)") + ggtitle(title)
}

npg <- 20
mu <- 2^c(8,10,9,11,10,12)
cond <- rep(rep(c("A","B"),each=npg),3)
geno <- rep(c("I","II","III"),each=2*npg)
counts <- rnbinom(6*npg, mu=rep(mu,each=npg), size=1/.01)
d <- data.frame(log2c=log2(counts+1), cond, geno)

Cairo(800, 800, "desvig000.png", bg="white")
plotit2(d, "Gene1") + ylim(7,13)
dev.off()
