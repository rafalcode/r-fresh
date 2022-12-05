#!/usr/bin/env Rscript
#from lmFit's help ?lmFit
library(limma)
library(Cairo)

# Simulate gene expression data for 100 probes and 6 microarrays
# Microarray are in two groups
# First two probes are differentially expressed in second group
# Std deviations vary between genes with prior df=4
sd <- 0.3*sqrt(4/rchisq(100,df=4))
y <- matrix(rnorm(100*6,sd=sd),100,6)
rownames(y) <- paste("Gene",1:100)
y[1:2,4:6] <- y[1:2,4:6] + 2
design <- cbind(Grp1=1,Grp2vs1=c(0,0,0,1,1,1))
options(digits=3)
stop("o!")
# Ordinary fit
fit <- lmFit(y,design)
ebfit <- eBayes(fit) # "moderation" invokes some Bayesianism via informed prior for taming purposes.
tt <- topTable(ebfit,coef=2) # BTW, you can't run topTable without the eBayes() step.
# dim(fit) # 100 2
# colnames(fit) # same as colnames for design.
# rownames(fit)[1:10] # gene names
# the fit object can have alot of components .. even head() will shoot off screen
# the many compoenents can be named via:
# names(fit) # the components of the object ebfit.

# Fold-change thresholding
fit2 <- treat(ebfit,lfc=0.1)
ttt <- topTreat(fit2,coef=2)

# Volcano plot
CairoPNG("lmf0.png", 800, 800)
volcanoplot(ebfit,coef=2,highlight=2)
dev.off()

# Mean-difference plot
CairoPNG("lmf1.png", 800, 800)
plotMD(ebfit,column=2)
dev.off()

# Q-Q plot of moderated t-statistics
CairoPNG("lmf2.png", 800, 800)
qqt(ebfit$t[,2], df = ebfit$df.residual + ebfit$df.prior) # qqt() a limma func.
abline(0,1)
dev.off()

# Various ways of writing results to file
write.fit(ebfit,file="lmf0.txt")
write.table(ebfit,file="lmf1.txt")

# Fit with correlated arrays
# Suppose each pair of arrays is a block
block <- c(1,1,2,2,3,3)
dupcor <- duplicateCorrelation(y,design,block=block)
# dupcor$consensus.correlation
fit3 <- lmFit(y, design, block=block, correlation=dupcor$consensus)
# don't know why they leave the following out:
# dim(fit3)
fit3 <- eBayes(fit3)
tt3 <- topTable(fit3,coef=2)

# Fit with duplicate probes # designed as more interesting? i.e. eBayes applied.
# Suppose two side-by-side duplicates of each gene
rownames(y) <- paste("Gene",rep(1:50,each=2))
dupcor2 <- duplicateCorrelation(y,design,ndups=2)
# dupcor2$consensus.correlation
fit4 <- lmFit(y,design,ndups=2,correlation=dupcor2$consensus)
# dim(fit4)
fit4 <- eBayes(fit4)
tt4 <- topTable(fit4,coef=2)
