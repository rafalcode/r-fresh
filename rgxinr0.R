#!/usr/bin/env Rscript
# A common question in R: how doe si tmanage regex?
# the answer is "badly"
# it will avoid saying whether there was not a match

cw <- getwd()
# gsub("(.+)/.+$", "\\1", cw, perl=T)

# OK the typical underscore parsing
n <- "20210204_TJC001_709"
p2 <- gsub("^[^_]+_(.+)$", "\\1", n, perl=T)

# a=can also eb achivee with strsplit
ss <- strsplit(n, "_")


p3 <- gsub("^\\d+_(.+)$", "\\1", n, perl=T)
if(p3 == "") {
    cat("p3 was null\n")
}
if(is.null(p3)) {
    cat("p3 was null\n")
}
n2 <- "20d10204_TJC001_709"
p4 <- gsub("^(\\d+)_(.+)$", "\\2", n2, fixed=T)
# p4 <- gsub("^(\\d+).+", "\\1", n2)
str(p4)
g <- grep("^\\d{4}", n)
str(g)

#A major encoutner .. the defautl for the non match is to give back the string as is!

nn <- c("20d10204_TJC001_709", "20d10204_TJC001_999")
p5 <- gsub("^[^_]+_(.+)$", "\\1", nn)
p6 <- gsub("^([^_]+)_.+", "\\1", nn)
p6
