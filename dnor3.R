# thi sis an important once because it has a histogram (made up count data)
# coexisting with a density plot
# otheise histogram is way too big.
library(ggplot2)
library(Cairo)
# Generate skewed data (not normal!)
data <- rexp(100, rate = 1)  # Exponential distribution

df <- data.frame(values = data)

CairoPNG("dnor3.png", 800, 800)
ggp <- ggplot(df, aes(x = values)) +
  geom_histogram(aes(y = after_stat(density)), 
                 bins = 30, 
                 fill = "lightblue", 
                 color = "black",
                 alpha = 0.7) +
  geom_density(fill = "red", alpha = 0.3, color = "red", linewidth = 1) +
  labs(title = "Histogram with Empirical Density Overlay") +
  theme_minimal()
show(ggp)
dev.off()
