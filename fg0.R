#!/usr/bin/env Rscript
# this script does what?
# fgsea .. time series RNA
library(Cairo)
library(GEOquery)
library(limma)

callgeo <- F
if(callgeo) {
    es <- getGEO("GSE200250", AnnotGPL = TRUE)[[1]]
} else {
    es <- readRDS("gse200250.rds")
}

# es, stands for expressionset, from Biobase
# it's a pretty big structure.
# es <- gse200250
es <- es[, grep("Th2_", es$title)]

es$time <- as.numeric(gsub(" hours", "", es$`time point:ch1`))

es <- es[, order(es$time)]

exprs(es) <- limma::normalizeBetweenArrays(log2(exprs(es)), method="quantile")

es <- es[order(rowMeans(exprs(es)), decreasing=TRUE), ]
es <- es[!duplicated(fData(es)$`Gene ID`), ]
rownames(es) <- fData(es)$`Gene ID`
es <- es[!grepl("///", rownames(es)), ]
es <- es[rownames(es) != "", ]

# fData, feature data
Biobase::fData(es) <- Biobase::fData(es)[, c("ID", "Gene ID", "Gene symbol")]

es <- es[head(order(rowMeans(exprs(es)), decreasing=TRUE), 12000), ]
# head(exprs(es))
# and so you get a 12000x18 "countset"
