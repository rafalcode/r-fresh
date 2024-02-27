#!/usr/bin/env Rscript
# this script does what? trying out methylmix
library(MethylMix)
library(Cairo)

# Load the three data sets needed for MethylMix
data(METcancer) # 14 rows (genes) and 251 TGCA.02 ... samples of some sort.
data(METnormal) # 14 genes x4 ... hmmm...? Well maybe only 4 normal samples!
data(GEcancer) #14 x 251

# Run methylmix on a small set of example data
res <- MethylMix(METcancer, GEcancer, METnormal)

# Plot the most famous methylated gene for glioblastoma
# i.e. you have to have identified it first.
CairoPNG("memix0.png", 800, 800)
mm <- MethylMix_PlotModel("MGMT", res, METcancer)
show(mm)
dev.off()

# Plot MGMT also with its normal methylation variation
CairoPNG("memix1.png", 800, 800)
MethylMix_PlotModel("MGMT", res, METcancer, METnormal = METnormal)
dev.off()

# Plot a MethylMix model for another gene
CairoPNG("memix2.png", 800, 800)
MethylMix_PlotModel("ZNF217", res, METcancer, METnormal = METnormal)
dev.off()

# Also plot the inverse correlation with gene expression (creates two separate plots)
CairoPNG("memix3.png", 800, 800)
MethylMix_PlotModel("MGMT", res, METcancer, GEcancer, METnormal)
dev.off()

# Plot all functional and differential genes
# well that could get quite voluminous
# for (gene in MethylMixResults$MethylationDrivers) {
#          MethylMix_PlotModel(gene, res, METcancer, METnormal = METnormal)
# }
