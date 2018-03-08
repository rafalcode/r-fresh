#!/usr/bin/env Rscript
# good exercise in splitting a string and then joining certain parts of it again.
library(stringi)
library(stringr)
cw <- getwd()
# str_match(cw, "(/([^/]+))+")
# cws<-str_split(cw, '/')
# nah, need the simplify
# cws<-str_split(cw, '/', simplify=T)
## we could the core strsplit, but it makes the paste render a vector, not clear how though.
cws<-strsplit(cw, '/')
# cwsp<-paste(cws, sep='/', collapse='')
# cwsp<-paste(cws[1:3], collapse='/')
cwsp<-paste(cws[1:3], sep='/')
# cwsp<-stri_join_list(cws[1:3], sep='/')
str(cws)
cws
str(cwsp)
cwsp

# Notes:
# stringr version
