#!/usr/bin/env Rscript
# checking out normalized coutns in DESeq2
# have to work out what factors actually means in sizefacor and normalizationfactors.
# are they weighting coeeficients of some kind.
library(Cairo)
library(DESeq2)
library(pheatmap)
library(RColorBrewer)

dds <- makeExampleDESeqDataSet(m=4)

# you get a deseq2 dataset obj from this
# you can see the count dat with assay(dds)
# and the pheno data with colData(dds)
# this could be a simulation of some sort, number of samples is 4 anyway, and two conditions.
dds <- estimateSizeFactors(dds)
# only now can we use the counts() function with "normalized" ON or T.
# the size factors will then appear as an extra column in the colData() func.
# let's keep conditions separate
# Conds <- data.frame(colData(dds)[,1], row.names=colnames(dds), col.names=colnames(colData(dds))[1])
Conds <- data.frame(condition=colData(dds)[,1], row.names=colnames(dds))
cou <-counts(dds, normalized=T)
# you get  a count matrix from this.

CairoPNG("denc0.png", 800, 800)
pheatmap(cou, show_rownames=F, cluster_cols=T, cluster_rows=F, scale="row",
#          color=colorRampPalette(rev(brewer.pal(n = 7, name ="RdBu"))), 
         color=colorRampPalette(rev(brewer.pal(n = 7, name ="RdBu")))(100), 
      clustering_distance_cols="euclidean", clustering_method="complete", border_color=FALSE,
      annotation_col=Conds)
# pheatmap(cou)
dev.off()

# don't use pheatmap's scale?
# bah comment out, because:
# yes you get almost the exact same, don't worry
# but probably better to use pheatmap's own.
# cous <- t(scale(t(cou)))
# CairoPNG("denc1.png", 800, 800)
# pheatmap(cous, show_rownames=F, cluster_cols=T, cluster_rows=F,
#       clustering_distance_cols="euclidean", clustering_method="complete", border_color=FALSE,
#       annotation_col=Conds)
# dev.off()
