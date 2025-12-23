#!/usr/bin/env Rscript
# this script does what?
library(ggplot2)
library(Cairo)

data <- rnorm(100, mean = 0, sd = 1)
df_data <- data.frame(values = data)

x <- seq(-4, 4, length.out = 200)
df_theory <- data.frame(x = x, y = dnorm(x, mean = 0, sd = 1))

CairoPNG("dnor2.png", 800, 800)
ggp <- ggplot() +
            geom_density(data = df_data, aes(x = values), 
                         fill = "lightblue", alpha = 0.5) +
            geom_line(data = df_theory, aes(x = x, y = y), 
                      color = "red", linewidth = 1) +
            # geom_histogram(data = df_data, stat="density", aes(x = values), bins=10, color="gray", alpha=0.25) +
            # above not working.
            labs(title = "Empirical vs Theoretical Density") +
            theme_minimal()
show(ggp)
dev.off()
