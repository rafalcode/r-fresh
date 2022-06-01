#!/usr/bin/env Rscript
library(Cairo)
library(gridExtra)
library(ggplot2)
library(ggrepel)
library(limma)

y <- matrix(rnorm(10*9),10,9)
colnames(y) <- paste0("s", 1:9)
pc0 <- prcomp(t(y))
xdf0 <- as.data.frame(pc0$x)

# perturb
y[,1:3] <- y[,1:3] + 5
# batch <- c("A","A","A","B","B","B","C","C","C")
batch <- c(rep("A", 3), rep("B", 6))
y2 <- removeBatchEffect(y, batch)

CairoPNG("boxplots.png", 1600, 800)
par(mfrow=c(1,2))
boxplot(as.data.frame(y),main="Original")
boxplot(as.data.frame(y2),main="Batch corrected")
dev.off()

pc1 <- prcomp(t(y))
xdf1 <- as.data.frame(pc1$x)
pc2 <- prcomp(t(y2))
xdf2 <- as.data.frame(pc2$x)
CairoPNG("pcaplots0.png", 1600, 800)
g0 <- ggplot(xdf0, aes(x=PC1, y=PC2, label=rownames(xdf0))) +geom_point() +geom_text_repel() +
    ggtitle("Original w/o batchfx") +
    theme(plot.title = element_text(color="purple", size=14, face="bold.italic"), axis.text=element_text(size=14), axis.title=element_text(size=14,face="bold"))
g1 <- ggplot(xdf2, aes(x=PC1, y=PC2, label=rownames(xdf2))) +geom_point() +geom_text_repel() +
    ggtitle("Batchfx applied and removed") +
    theme(plot.title = element_text(color="purple", size=14, face="bold.italic"), axis.text=element_text(size=14), axis.title=element_text(size=14,face="bold"))
grid.arrange(g0, g1, nrow=1)
dev.off()
CairoPNG("pcaplots2.png", 1600, 800)
g0 <- ggplot(xdf1, aes(x=PC1, y=PC2, label=rownames(xdf1))) +geom_point() +geom_text_repel() +
    ggtitle("Batchfx just applied") +
    theme(plot.title = element_text(color="purple", size=14, face="bold.italic"), axis.text=element_text(size=14), axis.title=element_text(size=14,face="bold"))
g1 <- ggplot(xdf2, aes(x=PC1, y=PC2, label=rownames(xdf2))) +geom_point() +geom_text_repel() +
    ggtitle("Batchfx applied and removed") +
    theme(plot.title = element_text(color="purple", size=14, face="bold.italic"), axis.text=element_text(size=14), axis.title=element_text(size=14,face="bold"))
grid.arrange(g0, g1, nrow=1)
dev.off()
