#!/usr/bin/env Rscript
# this is the very scattered MCSEA vignette
library(minfi)
library(FlowSorted.Blood.450k)
library(mCSEA)
library(Cairo)

data(mcseadata)
phenoTest$Pair <- rep(paste0("P", 1:10),2)
# myRank <- rankProbes(betaTest, phenoTest, refGroup = "Control")
myRank <- rankProbes(betaTest, phenoTest, refGroup = "Control", paired=T, pairColumn=3)
stop("ONLYME!")
myResults <- mCSEATest(myRank, betaTest, phenoTest, regionsTypes = "promoters", platform = "EPIC")
## Associating CpG sites to promoters
## Analysing promoters
## 38 DMRs found (padj < 0.05)
# mCSEATest() returns a list with the GSEA results and the association objects for
# each region type analyzed, in addition to the input data (methylation, phenotype
# and platform).
# ls(myResults)
## [1] "methData"
## [4] "promoters"

# "pheno"
# "platform"
# "promoters_association"
# 
# promoters is a data frame with the following columns (partially extracted from fgsea help):

#  pval: Estimated P-value.
#  padj: P-value adjusted by BH method.
#  log2err: Expected error for the standard deviation of the P-value logarithm.
#  ES: Enrichment score.
#  NES: Normalized enrichment score by number of CpGs associated to the
# eature.
#  size: Number of CpGs associated to the feature.
#  leadingEdge: Leading edge CpGs which drive the enrichment.

# head(myResults[["promoters"]][,-7])
##
# pval
# padj
# log2err
# ES
# NES size
## DMD
# 2.179762e-40 1.081162e-37 1.6533146 -0.9649723 -2.243737
# 65
## BANP
# 1.253191e-39 3.107915e-37 1.6407417 -0.9668041 -2.226523
# 59
# ## KTN1
# 2.382209e-16 1.476970e-14 1.0376962 -0.9605936 -1.964197
# 27
## XIAP
# 9.304952e-16 5.128062e-14 1.0175448 -0.9626065 -1.940590
## SEMA3B 4.231634e-15 1.908082e-13 0.9969862 -0.9603008 -1.935941
## GOLGB1 3.265203e-13 1.245800e-11 0.9325952 -0.9635271 -1.892163

# On the other hand, promoters_association is a list with the CpG probes associated to each feature:
# head(myResults[["promoters_association"]], 3)
## $YTHDF1
## [1] "cg18478105" "cg10605442" "cg27657131" "cg08514185" "cg13587582"
## [6] "cg25802399" "cg22485414" "cg03501095" "cg24092253" "cg12589387"
##
## $EIF2S3
## [1] "cg09835024" "cg06127902" "cg12275687" "cg00914804" "cg27345735"
## [6] "cg12590845" "cg25034591" "cg16712639" "cg07622257"
##
## $PKN3
## [1] "cg14361672" "cg06550760" "cg14204415" "cg11056832" "cg14036226"
## [6] "cg22365023" "cg20593100"

# You can also provide a custom association object between CpG probes and
# regions (customAnnotation option). This object should be a list with a structure
# similar to this:
# head(assocGenes450k, 3)
## $TSPY4
## [1] "cg00050873" "cg03443143" "cg04016144" "cg05544622" "cg09350919"
## [6] "cg15810474" "cg15935877" "cg17834650" "cg17837162" "cg25705492"
## [11] "cg00543493" "cg00903245" "cg01523029" "cg02606988" "cg02802508"
## [16] "cg03535417" "cg04958669" "cg08258654" "cg08635406" "cg10239257"
## [21] "cg13861458" "cg14005657" "cg25538674" "cg26475999"
##

## $TTTY14
## [1] "cg03244189" "cg05230942" "cg10811597" "cg13765957" "cg13845521"
## [6] "cg15281205" "cg26251715"
##
## $NLGN4Y
## [1] "cg03706273" "cg25518695" "cg01073572" "cg01498999" "cg02340092"
## [6] "cg03278611" "cg04419680" "cg05939513" "cg07795413" "cg08816194"
## [11] "cg09300505" "cg09748856" "cg09804407" "cg10990737" "cg18113731"
## [16] "cg19244032" "cg27214488" "cg27265812" "cg27443332"

# Step 3: Plotting the results
# Once you found some DMRs, you can make a plot with the genomic context of
# the interesting ones. For that, you must provide mCSEAPlot() function with
# the mCSEATest() results, and you must specify which type of region you want
# to plot and the name of the DMR to be plotted (e.g. gene name). There are
# some graphical parameters you can adjust (see mCSEAPlot() help). Take into
# account that this function connects to some online servers in order to get genomic
# information. For that reason, this function could take some minutes to finish
# the plot, specially the first time it is executed.
CairoPNG("mc0.png", 800, 800)
mCSEAPlot(myResults, regionType = "promoters", dmrName = "CLIC6", transcriptAnnotation = "symbol", makePDF = FALSE)
dev.off()
## Warning: Ensembl will soon enforce the use of https.
## Ensure the 'host' argument includes "https://"

# Chromosome 21
# 36.04 mb
# 
# 36.042 mb
# 36.041 mb
# 
# 36.043 mb
# 
# DNA Methylation
# 
# 0.8
# 0.6
# 0.4
# 0.2
# 
# Control
# 
# ENSEMBL annotation
# 
# KS
# leading
# edge
# 
# Case
# 
# You can also plot the GSEA results for a DMR with mCSEAPlotGSEA() function.
CairoPNG("mc1.png", 800, 800)
mCSEAPlotGSEA(myRank, myResults, regionType = "promoters", dmrName = "CLIC6")
dev.off()

# Integrating methylation and expression data
# If you have both methylation and expression data for the same samples, you can
# integrate them in order to discover significant associations between methylation
# changes in a DMR and an expression alterations in a close gene. mCSEAIntegrate() considers the DMRs identified by mCSEATest() passing a P-value
# threshold (0.05 by default). It calculates the mean methylation for each condition using the leading edge CpGs and performs a correlation test between
# this mean DMR methylation and the expression of close genes. This function
# automatically finds the genes located within a determined distance (1.5 kb) from
# the DMR. Only correlations passing thresholds (0.5 for correlation value and 0.05
# por P-value by default) are returned. For promoters, only negative correlations
# are returned due to this kind of relationship between promoters methylation
# and gene expression has been largely observed (Jones (2012)). On the contrary,
# only positive correlations between gene bodies methylation and gene expression
# are returned, due to this is a common relationship observed (Jones (2012)). For
# CpG islands and custom regions, both positive and negative correlations are
# returned, due to they can be located in both promoters and gene bodies.
# To test this function, we extracted a subset of 100 genes expression from bone
# marrows of 10 healthy and 10 leukemia patients (exprTest). Data was extracted
# from leukemiasEset package.
# Explore expression data
# head(exprTest, 3)

# Run mCSEAIntegrate function
resultsInt <- mCSEAIntegrate(myResults, exprTest, "promoters", "ENSEMBL")
## 0 genes removed from exprData due to Standar Deviation = 0
## Integrating promoters methylation with gene expression
# resultsInt
##

# Feature regionType
# 
# Gene Correlation
# 
# PValue
# 
# adjPValue
# 
# GATA2
# 
# promoters ENSG00000179348
# 
# -0.8908771 1.39373e-07 1.39373e-07
# 
# It is very important to specify the correct gene identifiers used in the expression
# data (geneIDs parameter). mCSEAIntegrate() automatically generates correlation plots for the significant results and save them in the directory specified by
# folder parameter (current directory by default).
# GATA2 promoter methylation vs. ENSG00000179348 expression
