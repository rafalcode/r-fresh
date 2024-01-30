#!/usr/bin/env Rscript
library(Cairo)

max.temp <- c(22, 27, 26, 24, 23, 26, 28)

CairoPNG("bp0.png", 800, 800)
barplot(max.temp)
dev.off()
