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
set.seed(42)
nrcou <- 4 # was: npg .. but I interpret: number of replicate counts for each condition-group
# mu <- 2^c(8,9,9,11,10,12) #rise due to condB in G3 is weaker.
mu <- c(256,512,512,2048,1024,1536)
# Note above mu, the power of two's style! These are six mu's which will help create 6 gene count vectors,
# one for each condition-group (2 conditions, 3 groups = 6)
phe <- data.frame(cond = rep(c("A","B"),nrcou*3),
                geno = rep(c("I","II","III"),each=2*nrcou))
rownames(phe) <- paste0("S", 1:(nrcou*6))

# here is the single gene to show DE
ndeg <- 3 # number of DE genes
# we'll need to modify the mu vector 
mu2 <- matrix(mu, byrow=T, ncol=2)
mu3 <- mu2 %x% rep(1, nrcou) # Kronecker to the rescue!
mu4 <- c(t(mu3))
counts <- matrix(rnbinom(6*nrcou*ndeg, mu=mu4, size=1/.01), nrow=ndeg, byrow=T)

# we'll want another #ng genes to show no DE but just small variations
ng <- 570000 # new genes
# OK different genes have very varied counts, these will be our seeds.
stg <- rnbinom(ng, mu=750, size=1)
# perturb the seed gene counts for condition A (untreated) 
e <- matrix(rnorm(ng*nrcou*3, sd=.05),nrow=ng) # different sample effect, no subgroup
e2 <- round(stg+stg*e) # perturbation is added
# for Cond B, just a bigger perturbation
f <- matrix(rnorm(ng*nrcou*3, sd=.2),nrow=ng) # condition effect
e3 <- round(e2+e2*f)
# note this approach will make some low coutns may be experience big changes, so that they 
# might occasionally even rival the true DEG's!

# Now we inteleave e and e3 because condA and CondB are interlevaed with each other.
# But R is very columns-first oriented ... especially c() so watch it here.
ee <-rbind(c(t(e2)), c(t(e3))) # fool c() to go row wise.
em <- matrix(c(ee), ncol=nrcou*6, byrow=T)
# now add the negligible effect genes to the DE gene.
# counts2 <- rbind(t(as.matrix(counts)), em)
counts2 <- rbind(counts, em)

logcou <- log2(counts2 + 1) # because lmFit() wants logs

# if you allow intercept (i.e. no ~0) you will get a logFC columns.
# if you don't, it appears it won't calculate difference for you.
des0 <- model.matrix(~cond*geno, data=phe)
fit0 <- lmFit(logcou, des0)
fit0 <- eBayes(fit0)
tt00 <- topTable(fit0)
tt01 <- topTable(fit0, coef=2)
tt02 <- topTable(fit0, coef=5)
tt03 <- topTable(fit0, coef=6)

# other designs
des <- model.matrix(~cond*geno, data=phe)
colnames(des) <- gsub(":", "_", colnames(des))
fit <- lmFit(logcou, des)

# Issue about Intercept being in bracked in the fit$coeeficients:
# colnames(fit$coefficients) <- gsub("[()]", "", colnames(colnames(fit$coefficients)))
cMat <- makeContrasts(condB_genoII - condB_genoIII, levels=des)
fit2 <- contrasts.fit(fit, cMat)
fit3 <- eBayes(fit2)
# fit30 <- eBayes(fit3)
tt3 <- topTable(fit3, coef=1)
# tt30 <- topTable(fit30)
