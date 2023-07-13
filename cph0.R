#!/usr/bin/env Rscript
# this script does what?
library(Cairo)
library(survival)
# canonicals from ?coxph

# Create the simplest test data set
test1 <- list(time=c(4,3,1,1,2,2,3),
              status=c(1,1,1,0,1,1,0),
              x=c(0,2,1,1,1,0,0),
              sex=c(0,0,0,0,1,1,1))

# Fit a stratified model
ct1 <- coxph(Surv(time, status) ~ x + strata(sex), test1)
# ct1 strucutre is a list of 20
# that firt arg is the formula .. i.e. Surv() has to be used, whereas lm() works with y ~ x
# strata() for one variable does nothing if it's a factor. with two, it paste0()'s the 2 of them

# Create a simple data set for a time-dependent model
test2 <- list(start=c(1,2,5,2,1,7,3,4,8,8),
              stop=c(2,3,6,7,8,9,9,9,14,17),
              event=c(1,1,1,1,1,1,1,0,0,0),
              x=c(1,0,0,1,0,1,1,1,0,0))
ct2 <- coxph(Surv(start, stop, event) ~ x, test2)
# summary(coxph(Surv(start, stop, event) ~ x, test2))

# Fit a stratified model, clustered on patients
bladder1 <- bladder[bladder$enum < 5, ]
cb <- coxph(Surv(stop, event) ~ (rx + size + number) * strata(enum), cluster = id, bladder1)

# Fit a time transform model using current age
cpage <- coxph(Surv(time, status) ~ ph.ecog + tt(age), data=lung, tt=function(x,t,...) pspline(x + t/365.25))
# pspline: yup splines are used.


# Cairo image template
# CairoPNG("fname.png", 800, 800)
# put plot command here
# dev.off()
