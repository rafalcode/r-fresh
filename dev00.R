#!/usr/bin/env Rscript
# from, hthe Interaction section of the DESEq2 vignette
# but echo=F for the RMD so this code doesn't come out.
# the "de" therefore means DESeq2

# Description we prepare a 6 valued vector called mu which will be the mean
# which is a parameter to the the rnbinom() and will generate count values around it.
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

# number of genes/rows.
# I got confused with this because of the long format.o# but it's 6 samples and 20 genes, which is reasonable.
nrcou <- 20 # was: npg .. but I interpret: number of replicate counts for each condition-group
mu <- 2^c(8,10,9,11,10,12) # Note the power of two's style! These are six mu's which will help create 6 gene count vectors, one for each condition-group (2 conditions, 3 groups = 6)
cond <- rep(rep(c("A","B"),each=nrcou),3)
geno <- rep(c("I","II","III"),each=2*nrcou)
# table(cond, geno)
# we prepare a counts vector in logn form, i.e. not like the usual matrix
# which in this case would be 20x6. Long form is continuous, so 120 long, filling out each column 
counts <- rnbinom(6*nrcou, mu=rep(mu,each=nrcou), size=1/.01)
couma <- matrix(counts, ncol=6)
d <- data.frame(log2c=log2(counts+1), cond, geno)
l <- lm(log2c ~ cond, data=d)

# Now let's change things, 
mu[4] <- 2^12 # so group2 for 2nd condition is higher.
mu[6] <- 2^8
counts2 <- rnbinom(6*nrcou, mu=rep(mu,each=nrcou), size=1/.01) # new count set based on the modified mu's.
couma2 <- matrix(counts2, ncol=6)
d2 <- data.frame(log2c=log2(counts2 + 1), cond, geno)
l2 <- lm(log2c ~ cond, data=d2)
l22 <- lm(log2c ~ cond + geno + geno:cond, data=d2)
