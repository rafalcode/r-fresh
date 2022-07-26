#!/usr/bin/env Rscript
# Kasper is the minfo author and dev, much given to methylation, so it's interesting to see
# him here talking about limma, though minfi may indeed usually use 
# limma for later analysis. In any case. Hansen has collaborated with Irizarry in the past
# ref. https://kasperdanielhansen.github.io/genbioconductor/html/limma.html
library(limma)
library(leukemiasEset)
data(leukemiasEset)


# Quick Leuk primer: there are 4 main types: 
# 1. Acute Lymphoblastic Leukemia (ALL). Subtype: c-ALL / pre-B-ALL without t(9;22)
# 2. Acute Myeloid Leukemia (AML). Subtype: Normal karyotype
# 3. Chronic Lymphocytic Leukemia (CLL)
# 4. Chronic Myeloid Leukemia (CML)
the NoL
# Let us ask which genes are differentially expressed between the ALL type and normal controls. First we subset the data and clean it up
ourData <- leukemiasEset[, leukemiasEset$LeukemiaType %in% c("ALL", "NoL")]

# so interestingly when a contrast between 5 levels could be done, he opts for the easy one.
# Natural you could say, it's a tutorial
# but a 5 level comparison could be done.

ourData$LeukemiaType <- factor(ourData$LeukemiaType)

design <- model.matrix(~ ourData$LeukemiaType)
fit <- lmFit(ourData, design)
fit <- eBayes(fit)
tt <- topTable(fit)
