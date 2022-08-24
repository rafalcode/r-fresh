#!/usr/bin/env Rscript
# this script does what?
library(ggplot2)
library(Cairo)
library(powsimR)


data("CELseq2_Gene_UMI_Counts")
batch <- sapply(strsplit(colnames(CELseq2_Gene_UMI_Counts), "_"), "[[", 1)
Batches <- data.frame(Batch = batch, stringsAsFactors = FALSE,
                        row.names = colnames(CELseq2_Gene_UMI_Counts))
data("GeneLengths_mm10")

# estimation
estparam_gene <- estimateParam(countData = CELseq2_Gene_UMI_Counts,
                               readData = NULL,
                               batchData = Batches,
                               spikeData = NULL,
                               spikeInfo = NULL,
                               Lengths = GeneLengths, MeanFragLengths = NULL,
                                RNAseq = 'singlecell', Protocol = 'UMI',
                                Distribution = 'NB', Normalisation = "scran",
                                GeneFilter = 0.1, SampleFilter = 3,
                                sigma = 1.96, NCores = NULL, verbose = TRUE)

# plotting
# plotParam(estParamRes = estparam_gene, Annot = T)

# define log fold change
p.lfc <- function(x) sample(c(-1,1), size=x,replace=T)*rgamma(x, shape = 1, rate = 2)
# set up simulations
setupres <- Setup(ngenes = 10000, nsims = 25, p.DE = 0.05, pLFC = p.lfc,
        n1 = c(48, 96, 384, 800), n2 = c(48, 96, 384, 800),
        Thinning = NULL, LibSize = 'equal',
        estParamRes = estparam_gene, estSpikeRes = NULL,
        DropGenes = TRUE, setup.seed = 5299, verbose = TRUE)

# Cairo image template
# CairoPNG("fname.png", 800, 800)
# put plot command here
# dev.off()
