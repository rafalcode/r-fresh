#!/usr/bin/env Rscript
# this script does what? DGCA .. their "basic" tutorial
library(ggplot2)
library(Cairo)
library(DGCA, quietly = TRUE)
# ref. https://cran.r-project.org/web/packages/DGCA/vignettes/DGCA_basic.html

data(darmanis)
data(design_mat)

# Note that the design matrix is a standard design matrix as used in other packages (e.g., limma, DESEq, MatrixEQTL), and specifies the group indices to be extracted from the original columns.
# To run the full differential correlation analysis and extract all of the top differentially correlated pairs, run this:

ddcor_res = ddcorAll(inputMat = darmanis, design = design_mat,
  compare = c("oligodendrocyte", "neuron"),
  adjust = "none", nPerm = 0, nPairs = 100)

ddcor_res = ddcorAll(inputMat = darmanis, design = design_mat,
  compare = c("oligodendrocyte", "neuron"),
  adjust = "perm", nPerm = 5, splitSet = "RTN4")
# head(ddcor_res)

# head(ddcor_res)
# Cairo image template
# CairoPNG("fname.png", 800, 800)
# put plot command here
# dev.off()
