# RF rendering of ATpoints Rmd doc.
# title: "https://www.biostars.org/p/461026/"
# author: "ATpoint"
# output: html_document
# obtained via wget https://gist.githubusercontent.com/ATpoint/c0d0e285cc94fa5c2087d4eb96f083c9/raw/f91fac93142012f84f41f0e3d96cfb24e3bb1269/biostars_461026.Rmd
library(data.table)
library(Cairo)
library(dplyr)
library(edgeR)
library(ggplot2)
library(limma)

parseProperly <- function(x){
  x <- x[,-c(ncol(x)-1,ncol(x))]
  x
}

# Read data drectly from GEO
m1 <- "https://ftp.ncbi.nlm.nih.gov/geo/series/GSE124nnn/GSE124167/suppl/GSE124167_FINAL_master_list_of_gene_counts_MIN.V4.txt.gz"
m2 <- "https://ftp.ncbi.nlm.nih.gov/geo/series/GSE124nnn/GSE124167/suppl/GSE124167_FINAL_master_list_of_gene_counts_MIN.sense.Pico.R.txt.gz"
m3 <- "https://ftp.ncbi.nlm.nih.gov/geo/series/GSE124nnn/GSE124167/suppl/GSE124167_FINAL_master_list_of_gene_counts_MIN.sense.Truseq.txt.gz"

method1 <- parseProperly(fread(m1, data.table=FALSE))
method2 <- parseProperly(fread(m2, data.table=FALSE))
method3 <- parseProperly(fread(m3, data.table=FALSE))

# Combine into a single table
counts <- Reduce(function(x, y) full_join(x, y, by="id"), list(method1, method2, method3))

# Replace NAs by zero
counts[is.na(counts)] <- 0

##Move "ids" to rownames to have a ure count matrix:
rownames(counts) <- counts$id
counts$id <- NULL

# Setup edgeR
theme_set(theme_bw()) #? R or edge R?

# Make a DGEList and add metadata
y <- DGEList(counts = counts)

y$samples$kit       <- c(rep("Pico", 6), rep("Illumina", 6), rep("V4", 6))
y$samples$treatment <- rep(c(rep("treated", 3), rep("untreated", 3)), 3)

# Normalize and obtain logcounts for QC
y <- calcNormFactors(y)
logCPMs <- cpm(y, log = TRUE)

# PCA

# Calculate rowwise variance
rv <- apply(logCPMs, 1, var)

# Sort decreasingly and take top 1000
o <- order(rv, decreasing=TRUE)
top1000 <- head(o, 1000)

# From the logCPMs subset for the top-1000
logCPM_top1000 <- logCPMs[top1000,]

# Run PCA
pca <- prcomp(t(logCPM_top1000))

# Combine PCA coordinates with the metadata from the DGEList
to_plot <- data.frame(pca$x, y$samples)

# Calculate how many % of total variance is explained by each principal component
percentVar <- pca$sdev^2 / sum( pca$sdev^2 )*100

# We focus here on PC1 and PC2
use.pcs <- c(1,2)
labs <- paste0(paste0("PC", use.pcs, " - "), paste0("Var.expl = ", round(percentVar[use.pcs], 2), "%"))

ggplot(to_plot, aes(x=PC1, y=PC2, color=kit, shape=treatment)) + 
  geom_point(size=3) +
  xlab(labs[1]) + ylab(labs[2])

# Batch correction
# correct for the batch which here is the "kit", using limma for this, a pretty famous function.
batch <- factor(y$samples$kit)
logCPMs_corrected <- removeBatchEffect(logCPMs, batch = batch)

# repeat PCA as before, using the same genes
logCPM_corrected_top1000 <- logCPMs_corrected[top1000,]

# Run PCA
pca <- prcomp(t(logCPM_corrected_top1000))

# Combine PCA coordinates with the metadata from the DGEList
to_plot <- data.frame(pca$x, y$samples)

# Calculate how many % of total variance is explained by each principal component
percentVar <- pca$sdev^2 / sum( pca$sdev^2 )*100

# We focus here on PC1 and PC2
use.pcs <- c(1,2)
labs <- paste0(paste0("PC", use.pcs, " - "), paste0("Var.expl = ", round(percentVar[use.pcs], 2), "%"))

CairoPNG("at0.png", 800, 800)
ggplot(to_plot, aes(x=PC1, y=PC2, color=kit, shape=treatment)) + 
  geom_point(size=3) +
  xlab(labs[1]) + ylab(labs[2])
dev.off()

# Differential analysis
# hmmm.. was this ever just plain subtraction?

# Design accounting for kit (=batch)
design <- model.matrix(~kit+treatment, y$samples)

# QLF workflow from edgeR
y <- estimateDisp(y, design)
fit <- glmQLFit(y, design)

# see head(design) => the fourth column is treatment which is what we want to test
fit  <- glmQLFTest(fit, coef = 4)

# get stats as a data.frame
tt <- data.frame(topTags(fit, n=Inf))

# Classify genes into significantly up and down
tt_modified <- tt %>% mutate(status=factor(case_when(logFC>0 & FDR<0.05 ~ "up", logFC<0 & FDR<0.05 ~ "down", TRUE ~ "not.signif"), levels=c("up", "not.signif", "down")))

# MA-plot
CairoPNG("at1.png", 800, 800)
ggplot(tt_modified, aes(x=logCPM, y=logFC, color=status)) +
  geom_point(size=1) +
  scale_color_manual(values=c("firebrick", "grey", "dodgerblue"))
dev.off()

# Volcano (logFC vs -log10(pvalue -- I prefer FDR))
CairoPNG("at2.png", 800, 800)
ggplot(tt_modified, aes(x=logFC, y=-log10(FDR), color=status)) +
  geom_point(size=1) +
  scale_color_manual(values=c("firebrick", "grey", "dodgerblue"))
dev.off()
