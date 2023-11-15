#!/usr/bin/env Rscript
# this script does what?
library(ggplot2)
library(Cairo)
library(ReactomePA)
library(DOSE)

data(geneList, package="DOSE")

de <- names(geneList)[abs(geneList) > 1.5] # de, diff exp I expect.
## [1] "4312"  "8318"  "10874" "55143" "55388" "991"
# i.e. ENTREZIDs!

x <- enrichPathway(gene=de, pvalueCutoff = 0.05, readable=TRUE)

y <- gsePathway(geneList, pvalueCutoff=.2, pAdjustMethod="BH", verbose=F)
# returns a gsea

CairoPNG("reactome.png", 800, 800)
vp <- ReactomePA::viewPathway("E2F mediated regulation of DNA replication", readable = TRUE, foldChange = geneList)
show(vp)
dev.off()

CairoPNG("reactomeerdp.png", 800, 800)
# x is enrichResult so this may work ...
dp <- dotplot(x, showCategory=30) + ggtitle("dotplot for Reactome")
show(dp)
dev.off()

CairoPNG("reactomegrdp.png", 800, 800)
# y is gseaResult so this may work ...
dp <- dotplot(y, showCategory=30) + ggtitle("dotplot for Reactome")
# only 1 pvalue? Odd. not very satisfying.
show(dp)
dev.off()

# conclusions
# dotplot on the enrichResults is much better
# reactome needs entrezids!
