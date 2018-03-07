#!/usr/bin/env Rscript

cw <- getwd()
gsub("(.+)/.+$", "$1", cw, perl=T)
cw
