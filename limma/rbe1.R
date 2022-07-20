#!/usr/bin/env Rscript
# simulation paired analysis
library(Cairo)
library(gridExtra)
library(ggplot2)
library(ggrepel)
library(limma)
library(RUVSeq)
library(sva)

# a 12 column matrices where columns are paired (see interlace effect below).
# first in pair having mean 0, the other mean 2.
# i.e. corresponding to 
# eoi <- factor(rep(c("P", "M"), 6)) # effect of interest
y00 <- cbind(matrix(rnorm(10*6, 8, 2),10,6), matrix(rnorm(10*6, 12, 2),10, 6))
y0 <- y00[,c(rbind(1:6, 7:12))] # interlace effect 1,7,2,8 etc.

phe <- as.data.frame(cbind(rep(c("C", "T"), 6), paste0("P", rep(1:6, each=2))))
rownames(phe) <- paste0("S", 1:12)
colnames(phe) <- c("Cond", "Pair")
colnames(y0) <- rownames(phe)
pc0 <- prcomp(t(y0))
xdf0 <- as.data.frame(pc0$x)

# perturb
# y[,1:3] <- y[,1:3] + 5
batch <- factor(rep(1:6, each=2)) # paired is batch effect
y <- y0+outer(rep(1,10) ,rep(0:5, each=2)*2) # "matrix-clever" way of cumulatively adding 2 to each pair of columns in y0
rownames(y) <- paste0("G", 1:10)

CairoPNG("boxplots0.png", 1600, 800)
par(mfrow=c(1,2))
boxplot(as.data.frame(y0),main="Original")
boxplot(as.data.frame(y),main="Pair-wise perturbed")
dev.off()

# Let's see what sva does with this
MM <- model.matrix(~0+Cond, data=phe)
# n.sv <- num.sv(y, mod=MM, method="leek")
# Error above: in eigen(wm) : infinite or missing values in 'x'
SV <- sva(y, mod=MM)
# will work on its own, and gives one SV, but it's just 1 for S1 and zero for all the rest.

y2 <- removeBatchEffect(y, batch)

CairoPNG("boxplots1.png", 1600, 800)
par(mfrow=c(1,2))
boxplot(as.data.frame(y0),main="Original")
boxplot(as.data.frame(y2),main="Batch corrected")
dev.off()

pc1 <- prcomp(t(y))
xdf1 <- as.data.frame(pc1$x)
pc2 <- prcomp(t(y2))
xdf2 <- as.data.frame(pc2$x)
CairoPNG("pcaplots0.png", 1600, 800)
g0 <- ggplot(xdf0, aes(x=PC1, y=PC2, label=rownames(xdf0), col=batch)) +geom_point() +geom_text_repel() +
    ggtitle("Original w/pairfx but w/o batchfx") +
    theme(plot.title = element_text(color="purple", size=14, face="bold.italic"), axis.text=element_text(size=14), axis.title=element_text(size=14,face="bold"), legend.position="none")
g1 <- ggplot(xdf2, aes(x=PC1, y=PC2, label=rownames(xdf2), col=batch)) +geom_point() +geom_text_repel() +
    ggtitle("Batchfx applied and then removed") +
    theme(plot.title = element_text(color="purple", size=14, face="bold.italic"), axis.text=element_text(size=14), axis.title=element_text(size=14,face="bold"), legend.position="none")
grid.arrange(g0, g1, nrow=1)
dev.off()
CairoPNG("pcaplots2.png", 1600, 800)
g0 <- ggplot(xdf1, aes(x=PC1, y=PC2, label=rownames(xdf1), col=batch)) +geom_point() +geom_text_repel() +
    ggtitle("W/pairfx and batchfx") +
    theme(plot.title = element_text(color="purple", size=14, face="bold.italic"), axis.text=element_text(size=14), axis.title=element_text(size=14,face="bold"), legend.position="none")
g1 <- ggplot(xdf2, aes(x=PC1, y=PC2, label=rownames(xdf2), col=batch)) +geom_point() +geom_text_repel() +
    ggtitle("Batchfx applied and removed") +
    theme(plot.title = element_text(color="purple", size=14, face="bold.italic"), axis.text=element_text(size=14), axis.title=element_text(size=14,face="bold"), legend.position="none")
grid.arrange(g0, g1, nrow=1)
dev.off()
