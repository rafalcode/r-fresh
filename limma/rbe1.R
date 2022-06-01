#!/usr/bin/env Rscript
# simulation paired analysis
library(Cairo)
library(gridExtra)
library(ggplot2)
library(ggrepel)
library(limma)
library(RUVSeq)

# a 12 column matrices where columns are paired, first in pair having mean 0, the other mean 2.
# i.e. corresponding to 
# eoi <- factor(rep(c("P", "M"), 6)) # effect of interest
y00 <- cbind(matrix(rnorm(10*6),10,6), matrix(rnorm(10*6, 1),10, 6))
y0 <- y00[,c(rbind(1:6, 7:12))]
rm(y00)
# note hte rbin

colnames(y0) <- paste0("S", 1:12)
pc0 <- prcomp(t(y0))
xdf0 <- as.data.frame(pc0$x)

# perturb
# y[,1:3] <- y[,1:3] + 5
batch <- factor(rep(1:6, each=2)) # paired is batch effect
y <- y0+outer(rep(1,10) ,rep(0:5, each=2)*8) # matrix -clever way of cumulatively adding 2 to each pair of columns in y0
CairoPNG("boxplots0.png", 1600, 800)
par(mfrow=c(1,2))
boxplot(as.data.frame(y0),main="Original")
boxplot(as.data.frame(y),main="Pair-wise perturbed")
dev.off()

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
