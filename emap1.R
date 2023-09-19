# emap0 was for a n EnrichRes object, this is for a compareClust
library(clusterProfiler)
library(enrichplot)
library(DOSE)
library(Cairo)
library(org.Hs.eg.db)

data(gcSample) # list of 8 genesets, all in entrez of course.

xx <- compareCluster(gcSample, fun="enrichGO", OrgDb="org.Hs.eg.db")

xx2 <- pairwise_termsim(xx) # sim refers to similarity
# you still get a compareLCuster result from this.

CairoPNG("emap1.png", 800, 800)
ema <- emapplot(xx2)
show(ema)
dev.off()
