#!/usr/bin/env Rscript
# https://rpubs.com/kylewbrown/r-colors
library(Cairo)

# 1.Define R Color Data ----
# RGB codes
color.rgb <- t(col2rgb(colors()))
# Hexadecimal codes
color.hex <- rgb(color.rgb[,1], color.rgb[,2], color.rgb[,3], maxColorValue = 255)
# Text highlighting
color.text <- ifelse(apply(color.rgb, 1, mean) > 127, "black", "white")
# Consolidate
color.df <- data.frame(name = colors(),
                       red = color.rgb[, "red"],
                       green = color.rgb[, "green"],
                       blue = color.rgb[, "blue"],
                       redp = color.rgb[, "red"]/255.,
                       greenp = color.rgb[, "green"]/255.,
                       bluep = color.rgb[, "blue"]/255.,
                       hex = color.hex,
                       text = color.text)

