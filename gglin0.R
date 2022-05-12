#!/usr/bin/env Rscript
# from https://r-charts.com/evolution/line-graph-multiple-lines-ggplot2/
library(ggplot2)
library(reshape)
library(Cairo)

# Consider the following data frame where each column represents the path of a brownian motion.
set.seed(2)

# Grid
t <- seq(0, 1, by = 0.01)
p <- length(t) - 1

# 5 paths
n <- 5
I <- matrix(rnorm(n * p, 0, 1 / sqrt(p)), n, p)

# Data frame
df1 <- data.frame(apply(I, 1, cumsum))

# In order to use your data frame in ggplot2 you will need to transform it into long format. You can achieve it making use of the melt function of the reshape package.
df <- data.frame(x = seq_along(df1[, 1]), df1)
# Long format
df <- melt(df, id.vars = "x")

# Given a data frame in long format like df it is possible to create a line chart with multiple lines in ggplot2 with geom_line the following way.

Cairo(800, 800, "ggl0.png", bg="white")
ggplot(df, aes(x = x, y = value, color = variable)) +geom_line()
dev.off()
 
# Lines width and style
# The styling of the lines can be changed making use of the arguments of geom_line, like linetype for changing the style of the line or lwd to change its width.
Cairo(800, 800, "ggl1.png", bg="white")
ggplot(df, aes(x = x, y = value, color = variable)) +geom_line(linetype = 3, lwd = 1.1)
dev.off()
 
# Color selection palette
cols <- c("#D43F3A", "#EEA236", "#5CB85C", "#46B8DA", "#9632B8")
Cairo(800, 800, "ggl2.png", bg="white")
ggplot(df, aes(x = x, y = value, color = variable)) +geom_line() +scale_color_manual(values = cols)
dev.off()
 
# Change the color of a line chart by group in ggplot2
# Note that using the previous method you can also highlight some lines of the chart, using the same color for all lines but some.

# Color selection
cols <- c("gray", "gray", "gray", "#5CB85C", "gray")

Cairo(800, 800, "ggl3.png", bg="white")
ggplot(df, aes(x = x, y = value, color = variable)) +geom_line() +scale_color_manual(values = cols)
dev.off()
 
# Change the title of the legend of a line chart in ggplot2
Cairo(800, 800, "ggl4.png", bg="white")
ggplot(df, aes(x = x, y = value, color = variable)) +geom_line() +guides(color = guide_legend(title = "My Title"))
dev.off()
 
# The labels of the legend can be modified making use of the labels argument of scale_color_discrete.
Cairo(800, 800, "ggl5.png", bg="white")
ggplot(df, aes(x = x, y = value, color = variable)) +geom_line() +scale_color_discrete(labels = paste("V", 1:5))
dev.off()
 
# Remove the legend of a line graph in ggplot2
# Just addYou can also get rid of the legend making use of legend.position = "none".
# +theme(legend.position = "none")
