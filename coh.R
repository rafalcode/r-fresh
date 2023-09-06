#!/usr/bin/env Rscript
# this script does what? An exercise with COHCAP (from the vignette)
library(ggplot2)
library(Cairo)
library(COHCAP)

dir = system.file("extdata", package="COHCAP")
beta.file = file.path(dir,"GSE42308_truncated.txt")
sample.file = file.path(dir,"sample_GSE42308.txt")

expression.file = file.path(dir,"expression-Average_by_Island_truncated.txt")
project.folder = tempdir()
cat("proj folder is", project.folder, "\n")
project.name = "450k_avg_by_island_test"

beta.table = COHCAP.annotate(beta.file, project.name, project.folder, platform="450k-UCSC")

filtered.sites = COHCAP.site(sample.file, beta.table, project.name, project.folder, ref="parental")

island.list = COHCAP.avg.by.island(sample.file, filtered.sites, beta.table,project.name, project.folder, ref="parental")

# apparently this next step integrates with gene expression:
COHCAP.integrate.avg.by.island(island.list, project.name, project.folder, expression.file, sample.file)
# and that's the last step in the vignette.

# Cairo image template
# CairoPNG("coh.png", 800, 800)
# put plot command here
# dev.off()
