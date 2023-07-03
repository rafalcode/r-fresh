#!/usr/bin/env Rscript
# this script does what?
library(ggplot2)
library(Cairo)
library(jtools)

states <- as.data.frame(state.x77)
fit1 <- lm(Income ~ Frost + Illiteracy + Murder + Population + Area + `Life Exp` + `HS Grad`, data = states, weights = runif(50, 0.1, 3))

fit2 <- lm(Income ~ Frost + Illiteracy + Murder + Population + Area + `Life Exp` + `HS Grad`, data = states, weights = runif(50, 0.1, 3))

fit3 <- lm(Income ~ Frost + Illiteracy + Murder + Population + Area + `Life Exp` + `HS Grad`, data = states, weights = runif(50, 0.1, 3))

CairoPNG("psumms.png", 800, 800)
# Plot all 3 regressions with custom predictor labels,
# standardized coefficients, and robust standard errors
ps <-plot_summs(fit1, fit2, fit3,
    # coefs = c("Frost Days" = "Frost", "% Illiterate" = "Illiteracy", "Murder Rate" = "Murder"),
    coefs = c("Frost", "Illiteracy", "Murder"),
    scale = TRUE, robust = TRUE)
show(ps)
dev.off()
