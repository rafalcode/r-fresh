#!/usr/bin/env Rscript
# this script does what?
library(ggplot2)
library(Cairo)

ourfile <- "de_df_for_volcano.rds"
ourpath <- "https://raw.githubusercontent.com/biocorecrg/CRG_RIntroduction/master/"
# download.file("https://raw.githubusercontent.com/biocorecrg/CRG_RIntroduction/master/de_df_for_volcano.rds", "de_df_for_volcano.rds", method="curl")
if(!file.exists(ourfile)) {
    download.file(paste0(ourpath, ourfile), ourfile, method="curl")
}
# Cairo image template
tmp <- readRDS("de_df_for_volcano.rds")

# remove rows that contain NA values
de <- tmp[complete.cases(tmp), ]

# CairoPNG("volc1.png", 800, 800)
# put plot command here
# dev.off()
