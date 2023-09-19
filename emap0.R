#!/usr/bin/env Rscript
# this script does what?
library(Cairo)
library(DOSE)

data(geneList)

de <- names(geneList)[1:100]

x <- enrichDO(de)
     x2 <- pairwise_termsim(x)
     emapplot(x2)
     # use `layout` to change the layout of map
     emapplot(x2, layout = "star")
   # use `showCategory` to select the displayed terms. It can be a number of a vector of terms.
     emapplot(x2, showCategory = 10)
     categorys <- c("pre-malignant neoplasm", "intestinal disease",
                    "breast ductal carcinoma")
     emapplot(x2, showCategory = categorys)

