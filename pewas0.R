#!/usr/bin/env Rscript
library(ggplot2)
library(Cairo)
library(pwrEWAS)
set.seed(1234)

# Running pwrEWAS by providing target maximal difference in methylation or by providing standard deviation of difference in methylation:
# rtdel = results_targetDelta
tT = "Blood adult"
# tarD: targetDelta
tarD <- c(0.02, 0.10, 0.15, 0.20)
rtdel <- pwrEWAS(minTotSampleSize = 20, maxTotSampleSize = 260, SampleSizeSteps = 40, NcntPer = 0.5,
                targetDelta = tarD, J = 100000, targetDmCpGs = 2500, tissueType = tT, 
                detectionLimit = 0.01, DMmethod = "limma", FDRcritVal = 0.05, core = 2, sims = 20)

# dSDvals = c(0.00390625, 0.02734375, 0.0390625, 0.052734375)
# results_deltaSD <- pwrEWAS(minTotSampleSize = 20, maxTotSampleSize = 260, SampleSizeSteps = 40, NcntPer = 0.5,
#                            deltaSD = dSDvals, J=100000, targetDmCpGs = 2500, tissueType = tT, detectionLimit = 0.01,
#                            DMmethod = "limma", FDRcritVal = 0.05, core = 2, sims = 20)
