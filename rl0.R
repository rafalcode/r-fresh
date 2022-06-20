#!/usr/bin/env Rscript
# this script does what?
library(ggplot2)
library(Cairo)
library(dplyr)
library(relayer)

df1 <- data.frame(x = c(0.5, 2.5, 4.5))

ggplot(df1, aes(x = x)) +
    geom_rect(aes(fill = x, y = 2)) +
    geom_tile(aes(fill2 = x, y = 1)) %>% rename_geom_aes(new_aes = c("fill" = "fill2")) +
    scale_colour_viridis_c(aesthetics = "fill", guide = "legend", name = "viridis A", option = "A") +
    scale_colour_viridis_c(aesthetics = "fill2", guide = "legend", name = "viridis D")
# Cairo image template
CairoPNG("fname.png", 800, 800)
# put plot command here
dev.off()
