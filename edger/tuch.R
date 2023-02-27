# tuch.R aka
# this is edgeR's example 4.1. in usersguide.
# actaully its matched pairs Oral Tumours (prob T) and Normal tissue (N).
library(edgeR)
library(org.Hs.eg.db)

rawdata <- read.delim("TableS1.txt", check.names=FALSE, stringsAsFactors=FALSE)
y <- DGEList(counts=rawdata[,4:9], genes=rawdata[,1:3])
# looks like 3 matched pairs, 15668 genes.

# 4.1.3 Annotation
# The study by Tuch et al. [40] was undertaken a few years ago, so not all of the RefSeq IDs
# provided by match RefSeq IDs currently in use. We retain only those transcripts with IDs in
# the current NCBI annotation, which is provided by the org.HS.eg.db package:
idfound <- y$genes$RefSeqID %in% mappedRkeys(org.Hs.egREFSEQ) # yep, this works.
y <- y[idfound,] # 15534 left.

# next bit is kind of a weird manoeuvre and not very suitable for a tutorial
# I suppose he wants to reverify that he's got hte right symbols.
# perhaps he didn't trust the SYMBOLs labels already there.
egREFSEQ <- toTable(org.Hs.egREFSEQ)
m <- match(y$genes$RefSeqID, egREFSEQ$accession)
y$genes$EntrezGene <- egREFSEQ$gene_id[m]
# Now use the Entrez Gene IDs to update the gene symbols:
egSYMBOL <- toTable(org.Hs.egSYMBOL)
m <- match(y$genes$EntrezGene, egSYMBOL$gene_id) # length(m) 15534!
y$genes$Symbol <- egSYMBOL$symbol[m]

# 4.1.4 Filtering and normalization
# Different RefSeq transcripts for the same gene symbol count predominantly the same reads.
# So we keep one transcript for each gene symbol. We choose the transcript with highest
# overall count:
o <- order(rowSums(y$counts), decreasing=TRUE)
y <- y[o,]
d <- duplicated(y$genes$Symbol) # returns logical vector.
y <- y[!d,] # this has major effect o the number of genes. Down to 10510.

# Why were there so many duplicates?
# well in fact the RefSeqIDs are differentA
# I used example of TTN and there are small differences only between them.
# However one looks very different with 312 exons instead of the usual 191/2, but actually one has 46 exons
# though that's quite a bizarre one (TTN is titin).
# Actually he does admit this:
# "Different RefSeq transcripts for the same gene symbol count predominantly the same reads."
# OK then he says that he will choose the one wiht highest count. OK, fair enough.

# Normally we would also filter lowly expressed genes (<50 reads). For this data however, all transcripts already
# have at least 50 reads for all samples of at least one of the tissues types.
# Recompute the library sizes:
y$samples$lib.size <- colSums(y$counts)

# Use Entrez Gene IDs as row names:
rownames(y$counts) <- rownames(y$genes) <- y$genes$EntrezGene
y$genes$EntrezGene <- NULL

# TMM normalization is applied to this dataset to account for compositional difference between the libraries.
y2 <- calcNormFactors(y)

# 4.1.5 Data exploration
# The first step of an analysis should be to examine the samples for outliers and for other
# relationships. The function plotMDS produces a plot in which distances between samples
# correspond to leading biological coefficient of variation (BCV) between those samples:

# plotMDS(y)

# 4.1.6 The design matrix
# Before we fit negative binomial GLMs, we need to define our design matrix based on the
# experimental design. Here we want to test for differential expression between tumour and
# normal tissues within patients, i.e. adjusting for differences between patients. In statistical
# terms, this is an additive linear model with patient as the blocking factor:
Patient <- factor(c(8,8,33,33,51,51))
Tissue <- factor(c("N","T","N","T","N","T"))
# data.frame(Sample=colnames(y),Patient,Tissue)
# Sample Patient Tissue

design <- model.matrix(~Patient+Tissue)
rownames(design) <- colnames(y2)

# This sort of additive model is appropriate for paired designs, or experiments with batch effects.

# 4.1.7 Estimating the dispersion
# We estimate the NB dispersion for the dataset.
y3 <- estimateDisp(y2, design, robust=TRUE)
# it's value is
# y$common.dispersion

# The square root of the common dispersion gives the coefficient of variation of biological
# variation. Here the common dispersion is found to be 0.159, so the coefficient of biological
# variation is around 0.4.
# The dispersion estimates can be viewed in a BCV plot:
# plotBCV(y)

# 4.1.8 Differential expression
# Now proceed to determine differentially expressed genes. Fit genewise glms:

fit <- glmFit(y3, design)

# Conduct likelihood ratio tests for tumour vs normal tissue differences and show the top genes:
lrt <- glmLRT(fit)
# topTags(lrt)

# Note that glmLRT has conducted a test for the last coefficient in the linear model, which we
# can see is the tumor vs normal tissue effect:
# colnames(design)

# The genewise tests are for tumor vs normal differential expression, adjusting for baseline
# differences between the three patients. The tests can be viewed as analogous to paired
# t-tests. The top DE tags have tiny p-values and FDR values, as well as large fold changes.
# Here’s a closer look at the counts-per-million in individual samples for the top genes:
o <- order(lrt$table$PValue)
#  cpm(y)[o[1:10],]

# We see that all the top genes have consistent tumour vs normal changes for the three patients.
# The total number of differentially expressed genes at 5% FDR is given by:
# summary(decideTests(lrt))

# Plot log-fold change against log-counts per million, with DE genes highlighted:
# plotMD(lrt)
# abline(h=c(-1, 1), col="blue")
# The blue lines indicate 2-fold changes.

# 4.1.9 Gene ontology analysis
# We perform a gene ontology analysis focusing on the ontology of biological process (BP).
# The genes up-regulated in the tumors tend to be associated with cell differentiation, cell
# migration and tissue morphogenesis:
go <- goana(lrt)
#  topGO(go, ont="BP", sort="Up", n=30, truncate=30)
# and that's the end right there.
