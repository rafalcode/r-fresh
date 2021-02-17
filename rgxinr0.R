#!/usr/bin/env Rscript
# A common question in R: how doe si tmanage regex?

cw <- getwd()
cw
gsub("(.+)/.+$", "\\1", cw, perl=T)
cw

# OK the typical underscore parsing
n <- "20210204_TJC001_709"
p2 <- gsub("^[^_]+_(.+)$", "\\1", n, perl=T)
p2


# a=can also eb achivee with strsplit
ss <- strsplit(n, "_")
ss
