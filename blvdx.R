#!/usr/bin/env Rscript
## Also from J Blischak .. this is the on about breast cancer
# but it's very sparingly documented and commented.
# ref. 
# https://jdblischak.github.io/dc-bioc-limma/vdx.html
library(ggplot2)
library(Cairo)
library(dplyr)
library(limma)
library(edgeR)
# Have to load Biobase after dplyr so that exprs function works
library(Biobase)
library(breastCancerVDX)
# i got his code off github
#
# Analysis of breast cancer VDX data for videos

# Prepare data -----------------------------------------------------------------

data("vdx")
# class(vdx)
# dim(vdx)
# head(pData(vdx))
# head(fData(vdx))
x <- exprs(vdx)
f <- fData(vdx)
p <- pData(vdx)

f <- f[, c("Gene.symbol", "EntrezGene.ID", "Chromosome.location")]
colnames(f) <- c("symbol", "entrez", "chrom")

# Recode er as 0 = negative and 1 = positive
p[, "er"] <- ifelse(p[, "er"] == 0, "negative", "positive")
p <- p[, c("id", "age", "er")]

# Explore data -----------------------------------------------------------------
CairoPNG("vdbox0.png", 800, 800)
boxplot(x[1, ] ~ p[, "er"], main = f[1, "symbol"])
dev.off()

eset <- ExpressionSet(assayData = x, phenoData = AnnotatedDataFrame(p), featureData = AnnotatedDataFrame(f))
# dim(eset)
CairoPNG("vdbox2.png", 800, 800)
boxplot(exprs(eset)[1, ] ~ pData(eset)[, "er"], main = fData(eset)[1, "symbol"])
dev.off()

# limma pipeline ---------------------------------------------------------------
design <- model.matrix(~er, data = pData(eset))
# head(design, 2)
# colSums(design)
# table(pData(eset)[, "er"])

fit <- lmFit(eset, design)
# head(fit$coefficients, 3)
fit <- eBayes(fit)
# head(fit$t, 3)
results <- decideTests(fit[, "erpositive"])
# summary(results)

# group-means ------------------------------------------------------------------
design <- model.matrix(~0 + er, data = pData(eset))
# head(design)
# colSums(design)

cm <- makeContrasts(status = erpositive - ernegative, levels = design)

fit <- lmFit(eset, design)
# head(fit$coefficients)
fit2 <- contrasts.fit(fit, contrasts = cm)
# head(fit2$coefficients)
fit2 <- eBayes(fit2)
results <- decideTests(fit2)
# summary(results)
# topTable(fit2)

# Visualize results  -----------------------------------------------------------
stats <- topTable(fit2, number = nrow(fit2), sort.by = "none")
# dim(stats)

# hist(runif(10000))

CairoPNG("vdhist.png", 800, 800)
hist(stats[, "P.Value"])
dev.off()

CairoPNG("vdvolc.png", 800, 800)
volcanoplot(fit2, highlight = 5, names = fit2$genes[, "symbol"])
dev.off()

# Enrichment -------------------------------------------------------------------

# topTable(fit2, number = 3)

# 1000 genes (10% in gene set), 100 are DE (10% in gene set)
fisher.test(matrix(c(10, 100, 90, 900), nrow = 2))

# 1000 genes (10% in gene set), 100 are DE (30% in gene set)
fisher.test(matrix(c(30, 100, 70, 900), nrow = 2))

head(fit2$genes, 3)
entrez <- fit2$genes[, "entrez"]

enrich_kegg <- kegga(fit2, geneid = entrez, species = "Hs")
topKEGG(enrich_kegg, number = 4)

enrich_go <- goana(fit2, geneid = entrez, species = "Hs")
topGO(enrich_go, ontology = "BP", number = 3)
