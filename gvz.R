#!/usr/bin/env Rscript
# this script does what? This is from Gviz vignette
# surprised I hadn't done this before
# did this at home but failed to igut push it.
# so gvi.R is almost the same.
library(Gviz)
library(Cairo)

# The most simple genomic features consist of start and stop coordinates, possibly overlapping each other. CpG islands or microarray probes are real life examples for this class of features. In the Bioconductor world those are most often represented as run-length encoded vectors, for instance in the IRanges and GRanges classes. To seamlessly integrate with other Bioconductor packages, we can use the same data structures to generate our track objects. A sample set of CpG island coordinates has been saved in the cpgIslands object and we can use that for our first annotation track object. The constructor function AnnotationTrack is a convenient helper to create the object.

library(GenomicRanges)

data(cpgIslands)
# class(cpgIslands)

# do you actually get 2 classes for thi?
# i.e.
## [1] "GRanges"
## attr(,"package")
## [1] "GenomicRanges"
# To clarify, the class name is "GRanges", and the "GenomicRanges" is just extra info
# saying which package GRanges comes from.
# in any case this is synthesised data, not real genome stuff. Correction, hg19 chr7 ... 
# In any case it's essnetialy 7 CpG islands (he says) but it's 10 really. right?

seqcg <- seqnames(cpgIslands) #an Rle class (run length encoding ... which comes from IRanges)
chr <- as.character(unique(seqcg))
gen <- genome(cpgIslands) # essential the Seqinfo slot, which is a class itself, so it's inside GRanges

atrack <- AnnotationTrack(cpgIslands, name = "CpG")
CairoPNG("gvz0.png", 800, 800)
plotTracks(atrack)
dev.off()

# so that was just the bars all you can guage is relative distance, scant info.
atrack <- AnnotationTrack(cpgIslands, name = "CpG")
CairoPNG("gvz0.png", 800, 800)
plotTracks(atrack)
dev.off()
