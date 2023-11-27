#!/usr/bin/env Rscript
# this script does what?
# have a go at directional with KEGG
library(Cairo)
library(clusterProfiler)
library(DOSE)

data(geneList, package="DOSE")
gene <- names(geneList)[abs(geneList) > 2]

kk <- enrichKEGG(gene         = gene,
                 organism     = 'hsa',
                 pvalueCutoff = 0.05)
