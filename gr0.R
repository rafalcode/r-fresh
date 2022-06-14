#!/usr/bin/env Rscript
# this script does what?
library(tidyverse)
library(Cairo)

n <- 10000
grad <- runif(n, min = 0, max = 100) %>% round()
x <- sample(letters, size = n, replace = T)

CairoPNG("fname.png", 800, 800)
# gg <- tibble(x, grad) %>% ggplot(aes(x = x, group = desc(grad), fill = grad)) + 
# gg <- tibble(x, grad) %>% ggplot(aes(x = x, group = grad, fill = grad)) + 
gg <- tibble(x, grad) %>% ggplot(aes(x = x, fill = grad)) + 
    geom_bar(stat = 'count') + 
    scale_fill_viridis_c()
show(gg)
dev.off()
