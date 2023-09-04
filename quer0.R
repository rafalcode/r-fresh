#!/usr/bin/env Rscript
# this script does what? Szymon Baluszek Chordoma study
library(GEOquery)

# gds0 <- getGEO("GDS507")

# gse0 <- getGEO("GSE781",GSEMatrix=F)

# GSE are series numbers, they're a container for
seriesnum <- "/mnt/sdb1/biodata/GSE230168/GSE230168" # Baluszek Chordoma study

retrieveit <- F
if(retrieveit) {
    gse <- getGEO("GSE230168", GSEMatrix=T)
    # actually T is the default. It's not about getting anything, it's about the parsing of the structures in the container.
    saveRDS(gse, paste0(seriesnum, "_geoqobj.rds"))
} else {
    gse <- readRDS(paste0(seriesnum, "_geoqobj.rds"))
}

# examine the object with show()

retrieveit <- F
if(retrieveit {
    gse2 <- getGEO("GSE230168", GSEMatrix=F)
    saveRDS(gse2, paste0(seriesnum, "_geoqobj2.rds"))
} else {
    gse2 <- readRDS(paste0(seriesnum, "_geoqobj2.rds"))

}

pd <- pData(phenoData(gse[[1]]))
# that will give you the meta data
# pd[,"methylation cluster:ch1"] gives the three found cluster

# GSMList doesn't work, gse turns out to be a list, of 1 actually. 1 expression set.

# sampleNames() will work o it. But they're GEO sample names.
# pd$title gives a bit unclear version otheir sample names
# pd$platform_id is the smae
# pd$geo_accession is the GSM name

# so I'm seeing that perhaps there is no expression?
# exprs(gse[[1]])
# will get you some sort of meth values. Yes, they go from 0 to 1, so it's mvalues.

# num [1:866238, 1:36] # not unuseful
