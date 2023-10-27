#!/usr/bin/env Rscript
# this script does what? DGCA .. their "basic" tutorial
library(ggplot2)
library(Cairo)
library(DGCA, quietly = TRUE)
# ref. https://cran.r-project.org/web/packages/DGCA/vignettes/DGCA_basic.html

data(darmanis)
# 'data.frame':   572 obs. of  158 variables
# those 572 rows are genenames
data(design_mat)
# each sample belongs to one of two categories neuron or oligodendrocyte.

# Note that the design matrix is a standard design matrix as used in other packages (e.g., limma, DESEq, MatrixEQTL), and specifies the group indices to be extracted from the original columns.
# To run the full differential correlation analysis and extract all of the top differentially correlated pairs, run this:

darmanis_filt = filterGenes(darmanis, filterTypes = c("central", "dispersion"), filterDispersionType = "cv", filterDispersionPercentile = 0.3) # 30th percentile, that is.

# for now continue tiwht unfiltered.
# getCors is jsut basic non-judgemental, and then we proceed to but we're more interested in differential, i.e. DCor
cor_res = DGCA::getCors(inputMat = darmanis, design = design_mat) 
dcPairs_res = DGCA::pairwiseDCor(cor_res, compare = c("oligodendrocyte", "neuron"))
# to get the adj
# "These are accessed during the extraction process of the table from the dcPairs class object"
# dd_pairs = dcTopPairs(dcPairs_res, nPairs="all", classify=T, adjust="BH")
dd_pairs = DGCA::dcTopPairs(dcPairs_res, nPairs=100, classify=T, adjust="BH")


# the filtered 
cor_res2 = DGCA::getCors(inputMat = darmanis_filt, design = design_mat) 
dcPairs_res2 = DGCA::pairwiseDCor(cor_res2, compare = c("oligodendrocyte", "neuron"))
# dd_pairs2 = DGCA::dcTopPairs(dcPairs_res2, nPairs="all", classify=T, adjust="BH")
# NO .. do not use "all" for nPairs here .. in other places perhaps, but not here.
dd_pairs2 = DGCA::dcTopPairs(dcPairs_res2, nPairs=37585, classify=T, adjust="BH")

