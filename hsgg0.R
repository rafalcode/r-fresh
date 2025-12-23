#!/usr/bin/env Rscript
# this script does what?
library(ggplot2)
library(Cairo)

x <- rnorm(100)
h <- hist(x, breaks=2, plot=F)

# muust create dataframe
df <- data.frame(
  xmin  = h$breaks[-length(h$breaks)],
  xmax  = h$breaks[-1],
  count = h$counts)

CairoPNG("hsgg0.png", 800, 800)
ggp <- ggplot(df) +
  geom_col(aes(x = xmin, y = count), width = df$xmax - df$xmin, align = "edge") +
  labs(x = "x", y = "Count") +
  theme_minimal()
show(ggp)
dev.off()
