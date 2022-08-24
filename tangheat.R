#!/usr/bin/env Rscript
# dave tang on pheatmap, says he likes it because of its annotating ability.
library(pheatmap)
library(DESeq2) # originally he had DESeq.
library(Cairo)
 
# if you can't install DESeq, I have hosted the file at https://davetang.org/file/TagSeqExample.tab
# example_file <- "https://davetang.org/file/TagSeqExample.tab"
 
# load data and subset
example_file <- "TagSeqExample.tab"
data <- read.delim(example_file, header=T, row.names="gene")
data_subset <- as.matrix(data[rowSums(data)>50000,])

CairoPNG("tang0.png", 800, 800)
pheatmap(data_subset, border=NA)
dev.off()
