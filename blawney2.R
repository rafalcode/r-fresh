#!/usr/bin/env Rscript
# insightfulpost!
# from https://support.bioconductor.org/p/129696/
library(DESeq2)
set.seed(0)

mock_dds <- makeExampleDESeqDataSet(n=10000, m=12) # m=12 is the default.
# without betaSD no DEGs.

# First extract the raw counts
raw_counts = counts(mock_dds)
# print('Original raw counts:')
# print(head(raw_counts))

# setup annotations
background = rep(c("sgCtrl","sgFoo"),each=6)
treatment = rep(rep(c("DMSO","T1"),each=3),2)
sample_annotations = cbind(background, treatment)
rownames(sample_annotations)<-paste0('Smp', 1:12)

# By setup, there should be no true effects so far.  
# For gene1, gene3, and gene5, add an effect due to the genetic background.
# To do this get the row mean for each of those genes and simply add that onto
# the sgFoo samples (for both DMSO and T1 treated conditions) 
selected_genes = c('gene1', 'gene3', 'gene5')
subset_means = floor(rowMeans(raw_counts[selected_genes,]))
# rowMeans output will give floats.

# create arrays for the samples we will eventually select:
sgFoo_samples = paste0('Smp', 7:12)
sgFoo_plus_treated_samples = paste0('Smp', 10:12)

# Add the perturbation to the sgFoo samples (sample7-sample12). Here, we add a 3-fold increase:
raw_counts[selected_genes, sgFoo_samples] = raw_counts[selected_genes, sgFoo_samples] + 3*subset_means
# has tje effect of raising each gene by 3 times its mean.

# Add yet another perturbation for only gene3 on the sgFoo+treated samples (sample10-sample12)
raw_counts['gene3', sgFoo_plus_treated_samples] = raw_counts['gene3', sgFoo_plus_treated_samples] + 3*subset_means['gene3']

# print('Perturbed raw counts:')
# print(head(raw_counts))

# Run DESEq2:
dds <- DESeqDataSetFromMatrix(countData = raw_counts,
                 colData = sample_annotations,
                 design = ~background+treatment+background:treatment)

dds <- DESeq(dds)
# print(resultsNames(dds))

# To get the background effect imposed on sgFoo--
res <- results(dds, name='background_sgFoo_vs_sgCtrl')
resOrdered <- res[order(res$pvalue),]
# print(head(resOrdered))
# 
# print('*****************')

# Check the interaction effect to see if gene3 pops up:
res <- results(dds, name="backgroundsgFoo.treatmentT1")
resOrdered <- res[order(res$pvalue),]
# print(head(resOrdered))
# print(resOrdered['gene3',])
