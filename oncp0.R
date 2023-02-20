#!/usr/bin/env Rscript
# oncoprint from the vignette
# ref. https://jokergoo.github.io/ComplexHeatmap-reference/book/oncoprint.html
library(ggplot2)
library(Cairo)
library(ComplexHeatmap)

mat <- read.table(system.file("extdata", package = "ComplexHeatmap", "tcga_lung_adenocarcinoma_provisional_ras_raf_mek_jnk_signalling.txt"), 
    header = TRUE, , sep = "\t") # the sep is necessary
# what we get in nmat is 172 obs. i.e. cases and 8 columns (variables) which are gene names.
# stop("o")
mat[is.na(mat)] <- "" # though there are plenty "empties" there are also "NA"
rownames(mat) <- mat[, 1] # this is because cases/obs are not in rownames yet, but in col1.
mat <- mat[, -1] # this gets rid of the first col, now that it's in rownames.
mat <-  mat[, -ncol(mat)]
# actually the above two are just quick ways to get rid of first and last cols.

mat <- t(as.matrix(mat))
# genes are rows now, and cases are columns .. that's the oncoprint style

col = c("HOMDEL" = "blue", "AMP" = "red", "MUT" = "#008000")
alter_fun = list(
    background = function(x, y, w, h) {
        grid.rect(x, y, w-unit(2, "pt"), h-unit(2, "pt"), 
            gp = gpar(fill = "#CCCCCC", col = NA))
    },
    # big blue
    HOMDEL = function(x, y, w, h) {
        grid.rect(x, y, w-unit(2, "pt"), h-unit(2, "pt"), 
            gp = gpar(fill = col["HOMDEL"], col = NA))
    },
    # big red
    AMP = function(x, y, w, h) {
        grid.rect(x, y, w-unit(2, "pt"), h-unit(2, "pt"), 
            gp = gpar(fill = col["AMP"], col = NA))
    },
    # small green
    MUT = function(x, y, w, h) {
        grid.rect(x, y, w-unit(2, "pt"), h*0.33, 
            gp = gpar(fill = col["MUT"], col = NA))
    }
)

# Cairo image template
# CairoPNG("fname.png", 800, 800)
# put plot command here
# dev.off()
