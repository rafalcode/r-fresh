#!/usr/bin/env Rscript
# this script does what?
library(ggplot2)
library(Cairo)
library(pheatmap)   
library(grid)

set.seed(123)

df<-data.frame( matrix(sample(30), ncol = 5))
colnames(df)<-LETTERS[1:5]
subj<-c("P1", "P2","P3", "T1", "T2","T3")
rownames(df)<-subj

aka2 = data.frame(ID = factor(rep(c("Pat","Trea"), each=3)))
rownames(aka2)<-subj
aka3 = list(ID = c(Pat = "white", Trea="blue"))

CairoPNG("ph2a.png", 800, 800)
pheatmap(t(scale(df)),
         annotation_col = aka2, 
         annotation_colors = aka3[1],
         annotation_legend = FALSE,
         gaps_col =  3,
         show_colnames = T, show_rownames = T, cluster_rows = F, 
         cluster_cols = F, legend = TRUE, 
         clustering_distance_rows = "euclidean", border_color = FALSE)

# Edit the relevant grob
# grid.ls(grid.force()) # "col_annotation" looks like it's the one to edit
# grid.gedit("col_annotation", gp = gpar(col="grey70"))
# grid.gedit("col_annotation", gp = gpar(col="black"))
dev.off()
