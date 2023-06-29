#!/usr/bin/env Rscript
# df() with no Cairo
library(jtools)

data(movies) # Telling R we want to use this data
fit <- lm(metascore ~ imdb_rating + log(us_gross) + genre5, data = movies)
su0 <- summ(fit)

fitg <- glm(metascore/100 ~ imdb_rating + log(us_gross) + genre5, data = movies,
            family = quasibinomial())

fit2 <- lm(metascore ~ imdb_rating + log(us_gross) + log(budget) + genre5, data = movies)

pdf("jt1.pdf", paper="a4r")
effect_plot(fitg, pred = imdb_rating, interval = TRUE, plot.points = TRUE, jitter = 0.05)
dev.off()

# pdf("jt2.pdf", paper="a4")
export_summs(fit, fit2, scale = TRUE, to.file="pdf", file.name="jt2.pdf")
# export_summs(fit, fit2, scale = TRUE, to.file="html", file.name="jt2.html")
# dev.off()
