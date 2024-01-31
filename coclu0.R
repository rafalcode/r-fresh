#!/usr/bin/env Rscript
# this script does what? using compareCLuster in clusterProfile package. there are three runs of it

# despite this being in the vignette, the plotting doesn't work.

# it was version 4.2.2
# he realises that  in v 4.8.2
library(Cairo)
library(clusterProfiler)

# NOTE:
# compareCluster() returns a comapreClusterResult object, from DOSE package
# a bit of a surprise, but DOSE must be within cluProf.

# unfort. there are errors in graphing.
# Error in as.double(y):
#   cannot coerce type 'S4' to vector of type 'double'
grapherrors <- F
# for v 4.2.2.

data(gcSample)
# this is a list of 8 lists. So, 8 gene lists or, better said, lists.
xx <- compareCluster(gcSample, fun="enrichKEGG", organism="hsa", pvalueCutoff=0.05)
#  as.data.frame(xx)

if(grapherrors) {
    CairoPNG("coclu0.png", 800, 800)
    plot(xx, type="dot", caption="KEGG Enrichment Comparison")
    dev.off()
}

## formula interface, another aspect
someentres <- c('1', '100', '1000', '100101467', '100127206', '100128071') # just some example entrez ids
mydf <- data.frame(Entrez=someentres, group = c('A', 'A', 'A', 'B', 'B', 'B'), othergroup = c('good', 'good', 'bad', 'bad', 'good', 'bad'))
xx.formula <- compareCluster(Entrez~group, data=mydf, fun='groupGO', OrgDb='org.Hs.eg.db')
# as.data.frame(xx.formula)

if(grapherrors) {
    CairoPNG("coclu1.png", 800, 800)
    plot(xx.formula, type="dot", caption="Entrez formula, GO Enrichment Comparison")
    dev.off()
}

## formula interface with more than one grouping variable
xx.formula.twogroups <- compareCluster(Entrez~group+othergroup, data=mydf, fun='groupGO', OrgDb='org.Hs.eg.db')
# as.data.frame(xx.formula.twogroups)

if(grapherrors) {
    CairoPNG("coclu2.png", 800, 800)
    plot(xx.formula.twogroups, type="dot", caption="Entrez formula 2 groups, GO Enrichment Comparison")
    dev.off()
}

# plot(xx, type="dot", caption="KEGG Enrichment Comparison")
CairoPNG("coclu3.png", 800, 800)
dp <- dotplot(xx) #xx unchanged.
show(dp)
dev.off()

## formula interface
mydf2 <- data.frame(Entrez=someentres, logFC = c(1.1, -0.5, 5, 2.5, -3, 3),
         group = c('A', 'A', 'A', 'B', 'B', 'B'),
         othergroup = c('good', 'good', 'bad', 'bad', 'good', 'bad'))
xx.formula <- compareCluster(Entrez~group, data=mydf2, fun='groupGO', OrgDb='org.Hs.eg.db')
# as.data.frame(xx.formula)
CairoPNG("coclu4.png", 800, 800)
dp <- dotplot(xx.formula)
show(dp)
dev.off()

xx.formula.twogroups <- compareCluster(Entrez~group+othergroup, data=mydf2, fun='groupGO', OrgDb='org.Hs.eg.db')
CairoPNG("coclu5.png", 800, 800)
dp <- dotplot(xx.formula.twogroups)
show(dp)
dev.off()

ck <- compareCluster(geneCluster = gcSample, fun = enrichKEGG)
CairoPNG("coclu6.png", 800, 800)
dp <- cnetplot(ck)
show(dp)
dev.off()
