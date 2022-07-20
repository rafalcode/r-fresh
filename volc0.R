#!/opt/R-4.1.3/bin/Rscript
# example from EnhancedVolcano's vignette.
# Cairo has no prolem generating the diagram.
library(Cairo)
library(airway)
library(magrittr)
library(org.Hs.eg.db)
library(DESeq2)
library(EnhancedVolcano)

data(airway)
airway$dex %<>% relevel('untrt')

# Annotate the Ensembl gene IDs to gene symbols:
ens <- rownames(airway)
symbols <- mapIds(org.Hs.eg.db, keys = ens, column = c('SYMBOL'), keytype = 'ENSEMBL')
symbols <- symbols[!is.na(symbols)] # everybody does this after ens to symbol changes. Those that don't convert, get deleted!
symbols <- symbols[match(rownames(airway), names(symbols))]
rownames(airway) <- symbols
keep <- !is.na(rownames(airway))
airway2 <- airway[keep,]

# Conduct differential expression using DESeq2 in order to create 2 sets of results:
dds <- DESeqDataSet(airway2, design = ~ cell + dex)
dds <- DESeq(dds, betaPrior=FALSE)
res <- results(dds, contrast = c('dex','trt','untrt'))
res <- lfcShrink(dds, contrast = c('dex','trt','untrt'), res=res, type = 'normal')

# 3.1 Plot the most basic volcano plot
# For the most basic volcano plot, only a single data-frame, data-matrix, or tibble of test results is required, containing point labels, log2FC, and adjusted or unadjusted P values. The default cut-off for log2FC is >|2|; the default cut-off for P value is 10e-6.
Cairo(800, 800, "vol0.png")
enh <- EnhancedVolcano(res, lab = rownames(res), x = 'log2FoldChange', y = 'pvalue')
show(enh)
dev.off()
