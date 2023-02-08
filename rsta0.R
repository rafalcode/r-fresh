#!/usr/bin/env Rscript
# this script does what?
# using rstan
# from https://www.r-bloggers.com/2019/01/an-introduction-to-stan-with-r/
# watch it though alot of things left unexplaine din this.
# despite fact it's one of the simplest prblems: Bernoulli
# seeing as they made an erro on the bar graph it's a degenerate tutorial.
library(rstan)
library(Cairo)

# STAN proposes you define a model first, out of program blocks, there can be six of these: data, transform_data, parameters, transformed_parameters, model, generated quantities.

# in this case, the call their modeldef bern.stan, ebcause this is a Bernouilli problem. 
# as you can see it's a string, just quite a ofrmally farranged string, that's all.
bern.stan <-
"data {
  int<lower=0> ntrials;               // number of trials
  int<lower=0, upper=1> y[ntrials];   // success on trial n
}

parameters {
  real<lower=0, upper=1> theta; // chance of success
}

model {
  theta ~ uniform(0, 1);        // prior
  y ~ bernoulli(theta);         // likelihood
}"

# Generate some simulated data
theta <- 0.3 # param for our simulation
ntrials <- 20 # lenght of data
y <- rbinom(ntrials, 1, theta) # ones and zeros, with .3 chance of a one.

# OK so despite that prob value, these 20 will have somethign different
# sum(y)/N would be how you estimate theta, which (they sya) you can call MLE. Maxlik estimator.
##  [1] 0 0 0 1 1 0 0 0 0 0 1 0 0 1 0 0 0 0 0 1
# in this one case it was 0.25, so fairly different to .3


# with all that done, we have the crowning (or crowding) operation, the "fit"
fit <- stan(model_code=bern.stan, data=list(y=y, N=ntrials), iter=5000)

# Now thwne that is finished, you can run
# print(fit, probs=c(0.1, 0.9))
#what's that probs? most likely it is a rstan specific overloading of print()
# of course the fit dstruct is pretty complicated, so you need to overload.

# they continue:
theta_draws <- extract(fit)$theta
# Calculating posterior mean (estimator)
# mean(theta_draws)
## [1] 0.2715866
# Calculating posterior intervals
# quantile(theta_draws, probs=c(0.10, 0.90))
##       10%       90% 
## 0.1569165 0.3934832
theta_draws_df <- data.frame(list(theta=theta_draws))
plotpostre <- ggplot(theta_draws_df, aes(x=theta)) +
  geom_histogram(bins=20, color="gray")

# Cairo image template
CairoPNG("pp.png", 800, 800)
show(plotpostre)
# put plot command here
dev.off()

# OK! in thre tutorial, bar chart centre on 0.5 which is way off!
# I'm sure it's a mistake.
# it should center around 0.3 and that's what I get.
