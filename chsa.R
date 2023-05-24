#!/usr/bin/env Rscript
# chsa .. complex heatmap single annotation
# however, this isn't one of those strip things.
# I thought that's what a single annotation was: a strop.
# obviously I'm wrong.
library(ComplexHeatmap)
library(Cairo)

mat = read.table(textConnection(
    "s1,s2,s3
    g1,snv;indel,snv,indel
    g2,,snv;indel,snv
    g3,snv,,indel;snv"), row.names = 1, header = TRUE, sep = ",", stringsAsFactors = FALSE)

mat = as.matrix(mat)
col = c(snv = "red", indel = "blue")

CairoPNG("chsa.png", 800, 800)
op <- oncoPrint(mat,
alter_fun = list(
    snv = alter_graphic("rect", width = 0.9, height = 0.9, fill = col["snv"]),
    indel = alter_graphic("rect", width = 0.9, height = 0.9, fill = col["indel"])
    ), col = col)
show(op)
dev.off()
