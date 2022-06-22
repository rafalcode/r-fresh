#!/usr/bin/env Rscript
# dmrcate's vignette
library(ExperimentHub)
library(missMethyl)
library(DMRcate)
library(Cairo)
library(IlluminaHumanMethylation450kanno.ilmn12.hg19)
library(IlluminaHumanMethylationEPICmanifest)

eh <- ExperimentHub()
FlowSorted.Blood.EPIC <- eh[["EH1136"]]
tcell <- FlowSorted.Blood.EPIC[,colData(FlowSorted.Blood.EPIC)$CD4T==100 |
colData(FlowSorted.Blood.EPIC)$CD8T==100]

detP <- detectionP(tcell)

remove <- apply(detP, 1, function (x) any(x > 0.01))
tcell <- preprocessFunnorm(tcell)

tcell <- tcell[!rownames(tcell) %in% names(which(remove)),]
tcellms <- getM(tcell)

nrow(tcellms)

tcellms.noSNPs <- rmSNPandCH(tcellms, dist=2, mafcut=0.05)

# Here we have 6 CD8+ T cell assays, and 7 CD4+ T cell assays; we want
tcell$Replicate[tcell$Replicate==""] <- tcell$Sample_Name[tcell$Replicate==""]
tcellms.noSNPs <- limma::avearrays(tcellms.noSNPs, tcell$Replicate)
tcell <- tcell[,!duplicated(tcell$Replicate)]
tcell <- tcell[rownames(tcellms.noSNPs),]
colnames(tcellms.noSNPs) <- colnames(tcell)
assays(tcell)[["M"]] <- tcellms.noSNPs
assays(tcell)[["Beta"]] <- ilogit2(tcellms.noSNPs)

# Next we want to annotate our matrix of M-values with relevant information.
type <- factor(tcell$CellType)
design <- model.matrix(~type)
myannotation <- cpg.annotate("array", tcell, arraytype = "EPIC", analysis.type="differential", design=design, coef=2)

# Now we can find our most differentially methylated regions with dmrcate().
dmrcoutput <- dmrcate(myannotation, lambda=1000, C=2)
results.ranges <- extractRanges(dmrcoutput, genome = "hg19")
groups <- c(CD8T="magenta", CD4T="forestgreen")
cols <- groups[as.character(type)]

Cairo(1000, 1000, "dm0plot.png", bg="white")
DMR.plot(ranges=results.ranges, dmr=1, CpGs=getBeta(tcell), what="Beta", arraytype = "EPIC", phen.col=cols, genome="hg19")
dev.off()

enrichment_GO <- goregion(results.ranges[1:100], all.cpg = rownames(tcell), collection = "GO", array.type = "EPIC")

enrichment_GO <- enrichment_GO[order(enrichment_GO$P.DE),]
