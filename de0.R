# from 
# http://bioconductor.org/packages/devel/bioc/vignettes/DESeq2/inst/doc/DESeq2.html#differential-expression-analysis
library(DESeq2)
library(pasilla)
library(Cairo)
library(ggplot2)

pasCts <- system.file("extdata", "pasilla_gene_counts.tsv", package="pasilla", mustWork=TRUE)
pasAnno <- system.file("extdata", "pasilla_sample_annotation.csv", package="pasilla", mustWork=TRUE)

cts <- as.matrix(read.csv(pasCts,sep="\t",row.names="gene_id"))
# cts[2,1] <- 140.48 # this will make DESDS() cough!
coldata <- read.csv(pasAnno, row.names=1)

coldata <- coldata[,c("condition","type")]
coldata$condition <- factor(coldata$condition)
coldata$type <- factor(coldata$type)
rownames(coldata) <- sub("fb", "", rownames(coldata))

# checks
# all(rownames(coldata) %in% colnames(cts))
## [1] TRUE
# all(rownames(coldata) == colnames(cts))
## [1] FALSE
cts <- cts[, rownames(coldata)]
# another check
# all(rownames(coldata) == colnames(cts))
## [1] TRUE
dds <- DESeqDataSetFromMatrix(countData = cts, colData = coldata, design = ~ condition)
dds <- DESeq(dds)

CairoPNG("dispests.png", 800, 800)
pp <- plotDispEsts(dds)
show(pp)
dev.off()
