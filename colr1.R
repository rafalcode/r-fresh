#!/usr/bin/env Rscript
# https://rpubs.com/kylewbrown/r-colors
# generate a c header file

# cat(paste0("typedef struct




# 1.Define R Color Data ----
# RGB codes
color.rgb <- t(col2rgb(colors()))
# Hexadecimal codes
color.hex <- rgb(color.rgb[,1], color.rgb[,2], color.rgb[,3], maxColorValue = 255)
# Text highlighting
color.text <- ifelse(apply(color.rgb, 1, mean) > 127, "black", "white")
# Consolidate
color.df <- data.frame(name = colors(),
                       redp = color.rgb[, "red"]/255.,
                       greenp = color.rgb[, "green"]/255.,
                       bluep = color.rgb[, "blue"]/255.)

for(i in 1:dim(color.df)[1]) {
    cat(paste0("{\"", color.df$name[i], "\", {", color.df$redp[i], ", ", color.df$greenp[i], ", ", color.df$bluep[i], "}}\n"))
}

