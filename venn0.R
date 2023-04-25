#!/usr/bin/env Rscript
# testing out VennDetail.
library(Cairo)
library(VennDetail)

# Hinder et al paper
data(T2DM)

# list to vennize
l2v <- list(Cortex = T2DM$Cortex$Entrez, SCN = T2DM$SCN$Entrez, Glom = T2DM$Glom$Entrez)
ven <- venndetail(l2v)

# Cairo image template
CairoPNG("venn0.png", 800, 800)
plot(ven)
dev.off()

# random two letter words:
ii <- 1+round(3*runif(64))
ss <- seq(1,64,2)
ss2 <- seq(2,64,2)
le <- paste0(letters[ii[ss]], letters[ii[ss2]])
l1 <- unique(le[1:16]) # it to unique-ize with sets.
l2 <- unique(le[17:32])

l2v <- list(Cortex = l1, SCN = l2)
ven <- venndetail(l2v)

# Cairo image template
CairoPNG("venn1.png", 800, 800)
plot(ven)
dev.off()
