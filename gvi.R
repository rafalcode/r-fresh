#!/usr/bin/env Rscript
# this script does what?
library(ggplot2)
library(Cairo)
library(Gviz)
library(GenomicRanges)

data(cpgIslands)
# class(cpgIslands)
## [1] "GRanges"
## attr(,"package")
## [1] "GenomicRanges"

chr <- as.character(unique(seqnames(cpgIslands)))
gen <- genome(cpgIslands) # simply "hg19"

atrack <- AnnotationTrack(cpgIslands, name = "CpG")

CairoPNG("gvi0.png", 800, 800)
plotTracks(atrack)
dev.off()

# note how the usual how() isn't necessary
# this has no ref ... all you're able to guage is the distance between bars
# only small value .. but necessary

# Well then let's try and attach that ref track
gtrack <- GenomeAxisTrack()
# Since a GenomeAxisTrack object is always relative to the other tracks that are plotted, there is little need for additional arguments. Essentially, the object just tells the plotTracks function to add a genomic axis to the plot. Nonetheless, it represent a separate annotation track just as the CpG island track does. We can pass this additional track on to plotTracks in the form of a list.

# ... and also the ideogram
itrack <- IdeogramTrack(genome = gen, chromosome = chr)

# Similar to the previous examples, we stick the additional track object into a list in order to plot it.

# CairoPNG("gvi1.png", 800, 800)
# plotTracks(list(itrack, gtrack, atrack))
# dev.off()
# it nows it's at around 46mbin chr7 because cpgIslands variable tells it.

data(geneModels)

grtrack <- GeneRegionTrack(geneModels, genome = gen, chromosome = chr, name = "Gene Model")

CairoPNG("gvi2.png", 800, 800)
plotTracks(list(itrack, gtrack, atrack, grtrack))
dev.off()

minbase <- grtrack@start
maxbase <- grtrack@end

# from https://stackoverflow.com/questions/70356629/how-to-increase-arrow-size-for-gvizs-ucsctrack-function-track
CairoPNG("gvi3.png", 800, 800)
rTrack <- UcscTrack(genome=gen, chromosome=chr, track="NCBI RefSeq", 
                    from=minbase, to=maxbase, trackType="GeneRegionTrack", 
                    rstarts="exonStarts", rends="exonEnds", gene="name", 
                    symbol="name2", transcript="name", strand="strand", 
                    fill="darkblue",stacking="squish", name="RefSeq", cex.title=2.8,
                    cex.axis=1.1, cex = 3, lwd = 2.5, showId=TRUE, geneSymbol=TRUE)

displayPars(rTrack) <- list(fontsize = 15) # changes title 
displayPars(rTrack) <- list(cex.group = 1.7) #  changes gene name 
#displayPars(rTrack) <- list(arrowHeadMaxWidth=30, arrowHeadWidth=30, shape="arrow") #  only changes the grey box
displayPars(rTrack) <- list(col.line = '#6e706b' ) # 
plotTracks(list(itrack, rTrack))
dev.off()
