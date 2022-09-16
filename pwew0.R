#!/usr/bin/env Rscript
# from pwrEWAS
# To execute the main pwrEWAS function the following to codes can be used. pwrEWAS
# allows the user to specify the effect size in one of two ways, by either providing a target
# maximal difference in methylation (“targetDelta”), or by providing the standard deviation
# of the simulated differnces (“deltaSD”). Only one of both arguments can be provided. If
# “targetDelta” is specified, pwrEWAS will automatically identify a standard deviation to simulate
# differences in methylation, such that the 99.99th percentile of the absolute value of simulated
# differences falls within a range around the targeted maximal difference in DNAm (see paper for
# additional details). If “deltaSD” is specified, pwrEWAS will simulate differences in methylation
# using the provided standard deviation (additional information provided in paper).
library("pwrEWAS")
# providing the targeted maximal difference in DNAm
res0 <- pwrEWAS(minTotSampleSize = 10,
                maxTotSampleSize = 50,
                SampleSizeSteps = 10,
                NcntPer = 0.5,
                targetDelta = c(0.2, 0.5),
                J = 100,
                targetDmCpGs = 10,
                tissueType = "Adult (PBMC)",
                detectionLimit = 0.01,
                DMmethod = "limma",
                FDRcritVal = 0.05,
                core = 4,
                sims = 50)

# providing the targeted maximal difference in DNAm
res2 <- pwrEWAS(minTotSampleSize = 10,
                maxTotSampleSize = 50,
                SampleSizeSteps = 10,
                NcntPer = 0.5,
                deltaSD = c(0.02, 0.05),
                J = 100,
                targetDmCpGs = 10,
                tissueType = "Adult (PBMC)",
                detectionLimit = 0.01,
                DMmethod = "limma",
                FDRcritVal = 0.05,
                core = 4,
                sims = 50)
