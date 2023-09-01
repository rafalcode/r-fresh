#!/usr/bin/env Rscript
library(DESeq2)

dds <- makeExampleDESeqDataSet(n=50000, m=36, betaSD=1) # betaSD=1 is more typical of my situation

# This is hard coded for only 2 conditions.
# they are called A and B, and the of the default 12 sample
# the first 6 are A and the second 6 are B. That simple.
# the model.matrix is also in there and it's condition_B_vs_A" 
# so differences "viewed from A" so to speak.
# the dds object already has condition built in .. with an intercept.

# big note:
# BetaSD control the magnitude of the diffexp between the two conditions.
# Like, you look at the function and none of that is explicit
dds <- DESeq(dds)
res0 <- results(dds)
res2 <- results(dds, lfcThreshold=1)
rdf0 <- as.data.frame(res0)
rdf2 <- as.data.frame(res2)

# I was playing around with this
# i.e. length(which(rdf0$padj<.05 & rdf0$log2FoldChange>2.))  etc.
# if you have decetn gene numbers , lfcThreshold will indded turn out to be conservative
# so you can guess why people keep it as a separate filter.

# actuall y with m=36 and n=20000 res2 seems to be actually much better.
# yep, verified with n=50k too.o
# some links on this issue here:
# https://support.bioconductor.org/p/101504/

# in fact two gurus on the subject say you really should put it in.
