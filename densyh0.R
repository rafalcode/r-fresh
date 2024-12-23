# Yann Holtz on density
# there is a density() function in  but he uses ggplot here.
# Libraries
library(ggplot2)
library(dplyr)
library(Cairo)

# Load dataset from github
# data <- read.table("https://raw.githubusercontent.com/holtzy/data_to_viz/master/Example_dataset/1_OneNum.csv", header=TRUE)
# Ehem, already problems this si a list of low precision fpnums, we can use readLines(). I also predownloaded.
# data <- readLines("1_OneNum.csv")
# well ok, yoiu need the header

# Make the histogram
CairoPNG("densyh0.png", 800, 800)
data %>%
  filter( price<300 ) %>%
  ggplot( aes(x=price)) +
    geom_density(fill="#69b3a2", color="#e9ecef", alpha=0.8)
dev.off()

