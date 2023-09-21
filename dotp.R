#!/usr/bin/env Rscript
# this script does what? The afamado dotplot
library(Cairo)
library(enrichplot)
library(DOSE)
library(clusterProfiler)
library(org.Hs.eg.db)
library(ggstar)

data(geneList)
de <- names(geneList)[1:100]

x <- enrichDO(de)
CairoPNG("dotp0.png", 800, 800)
dp <- dotplot(x)
show(dp)
dev.off()

# use `showCategory` to select the displayed terms. It can be a number of a vector of terms.
CairoPNG("dotp1.png", 800, 800)
dp <- dotplot(x, showCategory = 10)
show(dp)
dev.off()

#predefine the cats
categorys <- c("pre-malignant neoplasm", "intestinal disease", "breast ductal carcinoma", "non-small cell lung carcinoma")
CairoPNG("dotp2.png", 800, 800)
dp <- dotplot(x, showCategory = categorys)
show(dp)
dev.off()

# It can also graph compareClusterResult
data(gcSample)
data(gcSample)
xx <- compareCluster(gcSample, fun="enrichGO", OrgDb="org.Hs.eg.db")
xx2 <- pairwise_termsim(xx)

CairoPNG("dotp3.png", 800, 800)
dp <- dotplot(xx2)
show(dp)
dev.off()

CairoPNG("dotp4.png", 800, 800)
dp <- dotplot(xx2, shape = TRUE)
show(dp)
dev.off()

CairoPNG("dotp5.png", 800, 800)
dp <- dotplot(xx2, group = TRUE)
show(dp)
dev.off()

CairoPNG("dotp6.png", 800, 800)
dp <- dotplot(xx2, x = "GeneRatio", group = TRUE, size = "count")
show(dp)
dev.off()
