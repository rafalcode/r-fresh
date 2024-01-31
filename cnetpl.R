#!/usr/bin/env Rscript
# this script does what? using compareCLuster in clusterProfile package. there are three runs of it

# despite this being in the vignette, the plotting doesn't work.

# it was version 4.2.2
# he realises that  in v 4.8.2
library(Cairo)
library(clusterProfiler)
library(org.Hs.eg.db)

# NOTE:
# compareCluster() returns a comapreClusterResult object, from DOSE package
# a bit of a surprise, but DOSE must be within cluProf.

# unfort. there are errors in graphing.
# Error in as.double(y):
#   cannot coerce type 'S4' to vector of type 'double'
grapherrors <- F
# for v 4.2.2.

data(gcSample)

ck <- compareCluster(geneCluster = gcSample, fun = enrichKEGG)
CairoPNG("cnetpl0.png", 800, 800)
dp <- cnetplot(ck)
show(dp)
dev.off()

# note in his tutorial he has SYMBOL genenames come out
# but here it can  be seen that ENTREZIDs come out.

# ACTUALLY THE key is setReadable
ck <- setReadable(ck, OrgDb = org.Hs.eg.db, keyType="ENTREZID")
CairoPNG("cnetpl1.png", 800, 800)
dp <- cnetplot(ck)
show(dp)
dev.off()
