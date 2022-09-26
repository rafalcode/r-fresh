#!/usr/bin/env Rscript
# dave tang on pheatmap, says he likes it because of its annotating ability.
# note setting your own colors for the annotation bar is not easy.
library(pheatmap)
library(Cairo)
library(gplots)

# jet.colors <- colorRampPalette(c("#00007F", "blue", "#007FFF", "cyan", "#7FFF7F", "yellow", "#FF7F00", "red", "#7F0000"))
jet.colors <- colorRampPalette(c("#00007F", "#007FFF", "#7FFF7F", "#FF7F00", "#7F0000"))
# note special way you have to containerise the colours in
tisscol <- list(Tiss=c(Bx="Medium Orchid", BrM="Peru"))

# if you can't install DESeq, I have hosted the file at https://davetang.org/file/TagSeqExample.tab
# example_file <- "https://davetang.org/file/TagSeqExample.tab"
 
# load file with the data
dfile <- "sigbetas.csv"
data <- read.csv(dfile, row.names=1)
topn <- 500 # the top number to show (they're ordered by 
d <- data[1:topn,]
anfile <- "tatiss.csv"
ta <- read.csv(anfile, row.names=1)
colnames(ta) <- c("Tiss")
ta$Tiss <- factor(ta$Tiss)

CairoPNG("tang0.png", 800, 800)
pheatmap(d, main=paste0("Heatmap of top ", topn, " DMCpGs"), color=rich.colors(100),
         border=NA, show_rownames=F, annotation_col=ta, annotation_colors=tisscol,
         treeheight_row=0)
dev.off()
