#!/usr/bin/env Rscript
# from the DMRforPairs vignette ...
library(DMRforPairs)
data(DMRforPairs_data)

# head(CL.methy,2)
clm <- CL.methy
# this is a very handy dataframe. 711 obs and 14 variables.
# pretty much curated

# the first phase is tuning parameters, but you do have to inspect them
parameters <- expand.grid(min_distance = c(200,300), min_n = c(4,5))
# expand.gird is not in DMfor pkg, it's from S$Vectors, which masks R::base's version of this func.
resparams  <-  tune_parameters(parameters,
     classes_gene=clm$class.gene.related,
     classes_island=clm$class.island.related,
     targetID=clm$targetID,
     chr=clm$chromosome,
     position=clm$position,
     m.v=clm[,7:8],
     beta.v=clm[,11:12],
     recode=1,
     gs=clm$gene.symbol,
     do.parallel=2)

# these are all well known:
min_n=4 # min numb cpgs in region
d=200 # within which bp length?
dM=1.4 # minimum difference between the medians I think.
pval_th=0.05
experiment="DM4results_vignette"
method="fdr"
clr=c("red","blue","green")

# almost a repeat of the tuneparameters
output <- DMRforPairs(classes_gene=clm$class.gene.related,
                    classes_island=clm$class.island.related,
                    targetID=clm$targetID,
                    chr=clm$chromosome,
                    position=clm$position,
                    m.v=clm[,8:10],
                    beta.v=clm[,12:14],
                    min_n=4,
                    min_distance=200,
                    min_dM=1.4,
                    recode=1,
                    sep=";",
                    method="fdr",
                    debug.v=FALSE,
                    gs=clm$gene.symbol,
                    do.parallel=2)
