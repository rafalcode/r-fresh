#!/usr/bin/env Rscript
# good exercise in splitting a string and then joining certain parts of it again.
library(stringr)
cw <- getwd()
# str_match(cw, "(/([^/]+))+")
# cws<-str_split(cw, '/')
# cws<-str_split(cw, '/', simplify=TRUE)
# this here makes the paste render a vector
cws<-strsplit(cw, '/')
# cwsp<-paste(cws, sep='/', collapse='')
cwsp<-paste(cws[1:3], collapse='/')
cwsp
