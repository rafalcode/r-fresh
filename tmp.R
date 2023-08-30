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
