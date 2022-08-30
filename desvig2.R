#!/usr/bin/env Rscript
# http://bioconductor.org/packages/devel/bioc/vignettes/DESeq2/inst/doc/DESeq2.html
# This is an expansion, where extra non DE genes are simulated and added.
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

nrcou <- 5 # was: npg .. but I interpret: number of replicate counts for each condition-group
mu <- 2^c(8,10,9,11,10,11) #rise due to condB in G3 is weaker.
# Note above mu, the power of two's style! These are six mu's which will help create 6 gene count vectors,
# one for each condition-group (2 conditions, 3 groups = 6)
phe <- data.frame(cond = rep(c("A","B"),nrcou*3),
                geno = rep(c("I","II","III"),each=2*nrcou))
rownames(phe) <- paste0("S", 1:(nrcou*6))

# here is the single gene to show DE
ndeg <- 5 # number of DE genes
counts <- matrix(rnbinom(6*nrcou*ndeg, mu=rep(mu,nrcou), size=1/.01), nrow=ndeg, byrow=T)

# we'll want another #ng genes to show no DE but just small variations
ng <- 5000 # new genes
# OK different genes have very varied counts, these will be our seeds.
stg <- rnbinom(ng, mu=750, size=0.5)
# perturb the seed gene counts for condition A (untreated) 
e <- matrix(rnorm(ng*nrcou*3, sd=.05),nrow=ng) # different sample effect, no subgroup
e2 <- round(stg+stg*e) # perturbation is added
# for Cond B, just a bigger perturbation
f <- matrix(rnorm(ng*nrcou*3, sd=.2),nrow=ng) # condition effect
e3 <- round(e2+e2*f)

# Now we inteleave e and e3 because condA and CondB are interlevaed with each other.
# But R is very columns-first oriented ... especially c() so watch it here.
ee <-rbind(c(t(e2)), c(t(e3))) # fool c() to go row wise.
em <- matrix(c(ee), ncol=nrcou*6, byrow=T)
# now add the negligible effect genes to the DE gene.
# counts2 <- rbind(t(as.matrix(counts)), em)
counts2 <- rbind(counts, em)

logcou <- log2(counts2 + 1) # because lmFit() wants logs

des <- model.matrix(~0+cond*geno, data=phe)
colnames(des) <- gsub(":", "_", colnames(des))
fit <- lmFit(logcou, des)
cMat <- makeContrasts(condB-condA, levels=des)
fit2 <- contrasts.fit(fit, cMat)
fit3 <- eBayes(fit2)
fit30 <- eBayes(fit)
tt3 <- topTable(fit3)
tt30 <- topTable(fit30)
