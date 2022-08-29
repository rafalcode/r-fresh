#!/usr/bin/env Rscript
# generate beta distribution curves
# this is from David Robinson, he has a short article about beta distribution and baseball,
# and a longer article
library(Cairo)
library(ggplot2)
library(dplyr)

# df <- data.frame(a = c(81, 82, 81 + 100), b = c(219, 219, 219 + 200)) 
# first pair is mean, second lucky starter 4/10, third the veteran (300/1000)
df <- tibble(a = c(81, 81+4, 81 + 300), b = c(219, 219+10-4, 219+1000-300)) 

sim0 <- df %>% group_by(a, b) # this is a tibble function, so will convert df to a tibble.

sim <- sim0 %>%
       do(data_frame(x = seq(0, 1, .001), y = dbeta(x, .$a, .$b))) %>%
       mutate(Parameters = paste0("\u03B1 = ", a, ", \u03B2 = ", b)) %>%
       ungroup %>%
       mutate(Parameters = factor(Parameters, levels = unique(Parameters)))

CairoPNG("bb0.png", 800, 800)
sim %>% filter(a == 81) %>%
    ggplot(aes(x, y, color = Parameters)) + geom_line() +
    xlim(0, .5) + ylab("Density of beta") %>% show
dev.off()

CairoPNG("bb1.png", 800, 800)
sim %>% filter(a < 100) %>%
    ggplot(aes(x, y, color = Parameters)) + geom_line() +
    xlim(0, .5) + ylab("Density of beta") %>% show
dev.off()

CairoPNG("bb2.png", 800, 800)
sim %>% ggplot(aes(x, y, color = Parameters)) + geom_line() +
    xlim(0, .5) + ylab("Density of beta") %>% show
dev.off()
