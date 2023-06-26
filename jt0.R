#!/usr/bin/env Rscript
# this script does what? Experimenting with Jacob Long's jtools
# the main problem is the graphs
# how to incorporate Cairo into them?
library(Cairo)
library(jtools)

data(movies) # Telling R we want to use this data
fit <- lm(metascore ~ imdb_rating + log(us_gross) + genre5, data = movies)
su0 <- summ(fit)

fitg <- glm(metascore/100 ~ imdb_rating + log(us_gross) + genre5, data = movies,
            family = quasibinomial())

# Cairo image template
CairoPNG("jt0.png", 800, 800)
effect_plot(fitg, pred = imdb_rating, interval = TRUE, plot.points = TRUE, jitter = 0.05)
dev.off()
#CairoPDF("jt0.pdf", 0, 0, paper="a4r")
# CairoPDF("jt0.pdf", height=11.7, paper="a4r")
CairoPDF("jt0.pdf", paper="a4r")
effect_plot(fitg, pred = imdb_rating, interval = TRUE, plot.points = TRUE, jitter = 0.05)
dev.off()

CairoPNG("jt1.png", 800, 800)
plot_summs(fit)
dev.off()
