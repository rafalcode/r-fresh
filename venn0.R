#!/usr/bin/env Rscript
# testing out VennDetail.
library(VennDetail)
library(Cairo)
library(VennDetail)

# Hinder et al paper
data(T2DM)

ven <- venndetail(list(Cortex = T2DM$Cortex$Entrez, SCN = T2DM$SCN$Entrez, Glom = T2DM$Glom$Entrez))

# Cairo image template
CairoPNG("venn0.png", 800, 800)
plot(ven)
dev.off()
