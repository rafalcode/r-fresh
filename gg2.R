#!/usr/bin/env Rscript
# from https://r-charts.com/evolution/line-graph-multiple-lines-ggplot2/
# 5 random walks used as illurtation of geome_line and using repel labels.
library(ggplot2)
library(reshape)
library(Cairo)

# Consider the following data frame where each column represents the path of a brownian motion.
set.seed(2)

# Grid
t <- seq(0, 1, by = 0.01)
p <- length(t) - 1

# 5 random walks
n <- 5
I <- matrix(rnorm(n * p, 0, 1 / sqrt(p)), n, p)
# so I has 5 rows and 100 cols. i.e. 5 walks and 100 x-positions each.

# Data frame
df0 <- data.frame(apply(I, 1, cumsum)) # the usual random walk

# In order to use your data frame in ggplot2 you will need to transform it into long format. use reshape2's melt for this.
df1 <- data.frame(x = seq_along(df0[, 1]), df0)
df <- melt(df1, id.vars = "x")
# what do you get? well a df with three cols: "x", "variable" and "value"
# variable is X1 for all the values in column 1, followed (on the same column) by column2 with variable set to X2.

# Given a data frame in long format like df it is possible to create a line chart with multiple lines in ggplot2 with geom_line the following way.

Cairo(800, 800, "gg2.png", bg="white")
ggplot(df, aes(x = x, y = value, color = variable)) +geom_line()
dev.off()
