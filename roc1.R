#!/usr/bin/env Rscript
# this script does what? ROC curves from Revolutiona analytics, better I think
library(ggplot2)
library(Cairo)

set.seed(1)
sim_widget_data <- function(N, noise=100)
{
    # simulate good and bad widgets
  x <- runif(N, min=0, max=100)
  y <- 122 - x/2 + rnorm(N, sd=noise)
  bad_widget <- factor(y > 100)
  data.frame(x, y, bad_widget)
}
widget_data <- sim_widget_data(500, 10)

test_set_idx <- sample(1:nrow(widget_data), size=floor(nrow(widget_data)/4))

test_set <- widget_data[test_set_idx,]
training_set <- widget_data[-test_set_idx,]
library(ggplot2)
library(dplyr)
test_set %>% 
  ggplot(aes(x=x, y=y, col=bad_widget)) + 
  scale_color_manual(values=c("black", "red")) + 
  geom_point() + 
  ggtitle("Bad widgets related to x")
# Cairo image template
CairoPNG("fname.png", 800, 800)
# put plot command here
dev.off()
