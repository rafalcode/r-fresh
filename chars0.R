#!/usr/bin/env Rscript
# good exercise in splitting a string and then joining certain parts of it again.A
# so what different char/strings types does R actually have. That may be a simple question to answer
# but how do these same interaction with lists, vectors, matrices?
# What is a character matrix anyway?
library(stringi)
library(stringr)
cw <- getwd()
w2 <- c("tenuous", "tendentious")
cat("str(w2) is ")
str(w2)
cat("w2[2] is ")
w2[2]
cat("\n")

cws <- unlist(strsplit(cw, '/'))
cat("str(cws) is ")
str(cws)
cat("cws[2] is ")
cws[1:3]
cat("\n")
cwp <- paste(cws[1:3], collapse='/')
cwp
