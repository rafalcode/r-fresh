#!/usr/bin/env Rscript
# this script does what? tute of BEclear package
library(ggplot2)
library(Cairo)
library(BEclear)
library(knitr)

data("BEclearData")
# from this we get ex.data: a 250 x 40 matrix of betavals ... so there are 40 samples.
# adnd also ex.samples a 40x2 df giving sample_id and then 10 different batch_id
# so there is no discovery of unknown batches here.

batchEffect <- calcBatchEffects(
  data = ex.data, samples = ex.samples,
  adjusted = TRUE, method = "fdr")

# you get a list of 2 from this, and "med" which
# looks a bit weird, as each gene gets a different weighting for the same batch level.

#> INFO [2023-10-24 16:22:42] Transforming matrix to data.table
#> INFO [2023-10-24 16:22:42] Calculate the batch effects for 10 batches
#> INFO [2023-10-24 16:22:44] Adjusting p-values
mdifs <- batchEffect$med
pvals <- batchEffect$pval

summary <- calcSummary(medians = mdifs, pvalues = pvals)

score <- calcScore(ex.data, ex.samples, summary, dir = getwd())
#> INFO [2023-10-24 16:22:44] Generating a summary table
# knitr::kable(head(summary), caption = 'Summary over the batch affected gene-sample combination of the example data set')

# Cairo image template
# CairoPNG("fname.png", 800, 800)
# put plot command here
#dev.off()
