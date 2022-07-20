#!/usr/bin/env Rscript
# Gotta have heatmaps!
# hm3, alot of extraneous libs
library(ggplot2)
library(Cairo)
library(tidyverse)
library(hrbrthemes) # watch this guys wants extra fonts system wide!
library(viridis)
library(plotly)
# library(d3heatmap)
library(heatmaply)

# Load data 
# data <- read.table("https://raw.githubusercontent.com/holtzy/data_to_viz/master/Example_dataset/multivariate.csv", header=T, sep=";")
# colnames(data) <- gsub("\\.", " ", colnames(data))
# write.csv(data, "holtzy_multiv.csv", quote=F)
data <- read.csv("holtzy_multiv.csv", row.names=1)

# Select a few country
data <- data %>% 
  filter(Country %in% c("France", "Sweden", "Italy", "Spain", "England", "Portugal", "Greece", "Peru", "Chile", "Brazil", "Argentina", "Bolivia", "Venezuela", "Australia", "New Zealand", "Fiji", "China", "India", "Thailand", "Afghanistan", "Bangladesh", "United States of America", "Canada", "Burundi", "Angola", "Kenya", "Togo")) %>%
  arrange(Country) %>%
  mutate(Country = factor(Country, Country))

# Matrix format
mat <- data
rownames(mat) <- mat[,1]
mat <- mat %>% dplyr::select(-Country, -Group, -Continent)
mat <- as.matrix(mat)

# Heatmap
#d3heatmap(mat, scale="column", dendrogram = "none", width="800px", height="80Opx", colors = "Blues")

CairoPNG("fname.png", 800, 800)
p <- heatmaply(mat, 
        dendrogram = "none",
        xlab = "", ylab = "", 
        main = "",
        scale = "column",
        margins = c(60,100,40,20),
        grid_color = "white",
        grid_width = 0.00001,
        titleX = FALSE,
        hide_colorbar = TRUE,
        branches_lwd = 0.1,
        label_names = c("Country", "Feature:", "Value"),
        fontsize_row = 5, fontsize_col = 5,
        labCol = colnames(mat),
        labRow = rownames(mat),
        heatmap_layers = theme(axis.line=element_blank()))
dev.off()

