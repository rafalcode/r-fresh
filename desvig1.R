#!/usr/bin/env Rscript
# http://bioconductor.org/packages/devel/bioc/vignettes/DESeq2/inst/doc/DESeq2.html
# from, hthe Interaction section of the DESEq2 vignette
# but echo=F for the RMD so this code doesn't come out.

# Description we prepare a 6 valued vector called mu which will be the mean
# which is a parameter to the the rnbinom() and will generate count values around it.
library(ggplot2)
library(Cairo)
library(gridExtra)

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

nrcou <- 20 # was: npg .. but I interpret: number of replicate counts for each condition-group
mu <- 2^c(18,21,20,22,23,23)
# Note above mu, the power of two's style! These are six mu's which will help create 6 gene count vectors,
# one for each condition-group (2 conditions, 3 groups = 6)
cond <- rep(rep(c("A","B"),each=nrcou),3)
geno <- rep(c("I","II","III"),each=2*nrcou)
# table(cond, geno)
counts <- rnbinom(6*nrcou, mu=rep(mu,each=nrcou), size=1/.01)
d <- data.frame(log2c=log2(counts+1), cond, geno)

# a little experimentation with lm():
# lm() but for no intercept.
# l00 <- lm(log2c ~ 0+cond + geno, data=d)
# lm() with intercept
l0 <- lm(log2c ~ cond + geno, data=d)
l00 <- lm(log2c ~ cond, data=d)
# bove experiment, really shows there's very little difference between the two
# in fact the coefficents are just merely labelled differently.
# the first coeff is the first fact-level in l00 and "intercept" in l0.

# and back again to the experiment
l <- lm(log2c ~ cond + geno + geno:cond, data=d)

# Now let's change things, 
# mu[4] <- 2^12 # so group2 for 2nd condition is higher.
mu[6] <- 2^8
counts2 <- rnbinom(6*nrcou, mu=rep(mu,each=nrcou), size=1/.01) # new count set based on the modified mu's.
d2 <- data.frame(log2c=log2(counts2 + 1), cond, geno)
l2 <- lm(log2c ~ cond + geno + geno:cond, data=d2)
l20 <- lm(log2c ~ cond + geno, data=d2)
l200 <- lm(log2c ~ cond, data=d2)

# COmmentary up shot of the coefficients.
# BTW the coeeficients are clearly the most important part of lm(), they are the model.
# l gives
#  (Intercept)          condB         genoII        genoIII   condB:genoII  condB:genoIII
#       8.00484        2.01806        1.00432        1.92247       -0.07088        0.04398
# l2 gives
#   (Intercept)          condB         genoII        genoIII   condB:genoII  condB:genoIII
#        7.9635         2.0511         1.0464         2.0309         0.9261        -4.0792
# 
# So in l, you kind of have to forget about slope, these are categoircal not continuous variables
# the red lines suggest slope, but it's really just a higiher level of expression for the second level
# that's all. So condB is just 2 higher. higher than what? Well the implicit reference which is condA.
# the separate genoII and genoIII coeffs describe themselves relative to genoI and are independent of the conds.
# they either add or do not add.:w



# plots
# Cairo(1600, 800, "dev0_grp.png", bg="white")
# g1 <- plotit(d, "Gene 1") + ylim(7,13)
# g2 <- plotit(d2, "Gene 2") + ylim(7,13)
# grid.arrange(g1, g2, ncol=2)
# dev.off()

# what do the coefficients of lm's say when this is run?
#  (Intercept)          condB         genoII        genoIII   condB:genoII  condB:genoIII
#     8.01154        1.95778        1.01931        1.95538        0.07742        0.08589
#     8.0083         1.8564         0.9626         2.0095         1.2298        -3.8452
Cairo(1600, 800, "desvig0_nogrp.png", bg="white")
g1 <- plotit0(d, "Gene1") + ylim(7,13)
g2 <- plotit0(d2, "Gene2") + ylim(7,13)
grid.arrange(g1, g2, ncol=2)
dev.off()
