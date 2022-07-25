# ref. https://kasperdanielhansen.github.io/genbioconductor/html/limma.html
library(limma)
library(leukemiasEset)
data(leukemiasEset)

# Let us ask which genes are differentially expressed between the ALL type and normal controls. First we subset the data and clean it up
ourData <- leukemiasEset[, leukemiasEset$LeukemiaType %in% c("ALL", "NoL")]
ourData$LeukemiaType <- factor(ourData$LeukemiaType)

design <- model.matrix(~ ourData$LeukemiaType)
fit <- lmFit(ourData, design)
fit <- eBayes(fit)
tt <- topTable(fit)
