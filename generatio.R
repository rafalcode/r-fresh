#!/usr/bin/env Rscript
# this script does what?
# what is this gene ration thing in compareCluster()?
library(Cairo)
library(clusterProfiler)

genes <- letters[1:15]

gs_df <- data.frame("gs_name"=c(rep("genesetX", 10), rep("genesetY", 25)),
                    "entrez_gene"=c(letters[1:10], letters[2:26]))

enr <- clusterProfiler::enricher(gene = genes, TERM2GENE = gs_df, minGSSize=1)@result
