#!/usr/bin/env Rscript
# dave tang on pheatmap, says he likes it because of its annotating ability.
# note setting your own colors for the annotation bar is not easy.
library(ComplexHeatmap)
library(Cairo)
library(gplots)

# note special way you have to containerise the colours in (only for pheatmap)
# tisscol <- list(Tiss=c(Bx="Medium Orchid", BrM="Peru"))
tisscol <- c(Bx="Medium Orchid", BrM="Peru")
# paircol <- list(Pair=c(ERP="Dodger Blue", HER2P="Khaki", TNBC="Tomato"))
subcol <- c(ERP="Dodger Blue", HER2P="Khaki", TNBC="Tomato")
paircol <- c(P1="Steel Blue", P2="Rosy Brown", P3="Dark Sea Green", P4="Moccasin", P5="Plum", P6="Corn Silk", P7="Dark Salmon", P8="Spring Green", P9="Gold", P10="Powder Blue", P11="Maroon", P12="Lemon Chiffon", P13="Light Green", P14="Alice Blue", P15="Light Pink", P16="Sienna", P17="Aquamarine", P18="Gainsboro", P19="Dark Green", P20="Steel Blue", P21="Dark Orchid") 

# load file with the data
dfile <- "sigbetas.csv"
data <- read.csv(dfile, row.names=1)
topn <- 5000 # the top number to show (they're ordered by p-value)
da <- data[1:topn,]
dam <- as.matrix(da)

split = rep(1:4, times=c(8,8,7,7))

anfile <- "tatiss.csv"
ta <- read.csv(anfile, row.names=1)

# you have several annotatin bars, just append them with a comma,
# However the colours must be dealt with together.
ha <- HeatmapAnnotation(Tiss=ta$Tiss, Subprim=ta$Subprim, Pair=ta$Pair, border=T,
                        col=list(Tiss=tisscol, Subprim=subcol, Pair=paircol))
# ha0 <- columnAnnotation(Tiss=ta$Tiss, Subprim=ta$Subprim, Pair=ta$Pair)
CairoPNG("tangch.png", 800, 800)
hm <- Heatmap(dam, col=rich.colors(100), column_title=paste0("Heatmap of top ", topn, " (by p-value) DMCpGs"), 
         show_row_dend=F, show_row_names=F, bottom_annotation=ha, show_heatmap_legend=T,
         name="\U03B2-value", column_title_gp=gpar(fontface="bold",fontsize=18),
         column_split=split)
draw(hm, annotation_legend_side = "left")
dev.off()

# you can have the original ordering of the samples based on pairs being beside each other with:
# column_order=rownames(ta))
# however, the clustering will put similar overall intensities close to each other
# and columns that are far away as just more different to each other.
