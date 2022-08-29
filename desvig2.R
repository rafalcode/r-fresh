#!/usr/bin/env Rscript
# http://bioconductor.org/packages/devel/bioc/vignettes/DESeq2/inst/doc/DESeq2.html
# I wanted to expand the desvig1.R include more genes
library(ggplot2)
library(Cairo)
library(gridExtra)
library(limma)

plotit2 <- function(d, title) {
  # ggplot(d, aes(x=cond, y=log2c, group=1)) + geom_point() +
  ggplot(d, aes(x=cond, y=log2c)) + geom_point() +
    stat_summary(fun=mean, geom="line", colour="red", size=0.8) + 
    xlab("condition") + ylab("log2(counts+1)") + ggtitle(title)
}

plotit0 <- function(d, title) {
    # this is the orginal graphing func.
    # where the samples are jittered for easy visualisation
    # use the plotit2() function if you want an plainer less adorned graph.
  ggplot(d, aes(x=cond, y=log2c, group=geno)) + 
    geom_jitter(size=1.5, position = position_jitter(width=.15)) +
    facet_wrap(~ geno) + 
    stat_summary(fun=mean, geom="line", colour="red", size=0.8) + 
    xlab("condition") + ylab("log2(counts+1)") + ggtitle(title)
}

nrcou <- 6 # was: npg .. but I interpret: number of replicate counts for each condition-group
mu <- 2^c(8,10,9,11,10,11) #rise due to condB in G3 is weaker.
# Note above mu, the power of two's style! These are six mu's which will help create 6 gene count vectors,
# one for each condition-group (2 conditions, 3 groups = 6)
phe <- data.frame(sname = paste0("S", 1:(nrcou*6)),
                cond = rep(c("A","B"),nrcou*3),
                geno = rep(c("I","II","III"),each=2*nrcou))

counts <- rnbinom(6*nrcou, mu=rep(mu,nrcou), size=1/.01) # new count set based on the modified mu's.
d2 <- log2(counts + 1)

des <- model.matrix(~0+cond, data=phe)
des <- model.matrix(~0+cond, data=phe)
fit <- lmFit(counts, des)
cMat <- makeContrasts(condB-condA, levels=des)
fit <- contrasts.fit(fit, cMat)
fit <- eBayes(fit)
tt <- topTable(fit)
