#!/usr/bin/env Rscript
# do you know how to use dnorm?
# it gives you the density value of normal variate.
# by this I mean, is that if your variate was the mean itself (which actually is a variate)
# then you would get 0.39 or something which is the peak value of your distribution

# dnorm() always sticks to the theoretical, there is no variability in it 
# so it's useful for imposing the theoretical distri over you data.

library(ggplot2)
library(Cairo)

# You generated some data with rnorm()
data <- rbeta(500, 5, 2)
h <- hist(data)

# ggplot works best with data frames so we convert ...
df <- data.frame(values = data)

# Create a data frame for the thoeretical normal curve, first we set the x-range using seq
x_range <- seq(min(data), max(data), length.out = 20)
y_beta <- dbeta(x_range, 5, 2)
theo <- data.frame(x = x_range, y = y_beta)

# modify beta to reflect age categories rather than probabilities:
df <- df*100
theo$x <- theo$x*100

CairoPNG("dbe1.png", 800, 800)
ggp <- ggplot() +
       geom_col()
       geom_line(data = theo, aes(x = x, y = y), color = "red", linewidth = 1) +
       labs(title = "Histogram with Theoretical Beta Distribution", x = "Values", y = "Density") +
       theme_minimal()
show(ggp)
dev.off()

# this could have also have been done using ggplot's stat_function
# which seems to already incorporate the x range.
# stat_function(fun = dnorm,
#                args = list(mean = mean(data), sd = sd(data)),
#                color = "red",
#                linewidth = 1) +
