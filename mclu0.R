# the faithfuldataset, two vars, waiting and eruptions.
# to be clear the variables are for the durations, i.e. eruption duraction and waiting duration
# (so it's not really a poisson exercise, despite the word waitin.
# primary used to see if there's realtionship between er
library(Cairo)
set.seed(42)  # For reproducibility

n <- 272  # Same size as faithful dataset
# simerups, simualtion eruptions variable

simerups <- numeric(n)

# Estimate parameters:
short_mean <- mean(faithful$eruptions[faithful$eruptions < 3])
short_sd <- sd(faithful$eruptions[faithful$eruptions < 3])

long_mean <- mean(faithful$eruptions[faithful$eruptions >= 3])
long_sd <- sd(faithful$eruptions[faithful$eruptions >= 3])

p_short <- sum(faithful$eruptions < 3) / length(faithful$eruptions)  # Probability of short eruption
p_long <- 1 - p_short  # Probability of long eruption


# OK then let's go
for (i in 1:n) {
  if (runif(1) < p_short) {
    simerups[i] <- rnorm(1, mean=short_mean, sd=short_sd)
  } else {
    simerups[i] <- rnorm(1, mean=long_mean, sd=long_sd)
  }
}

# Plot simulated data
CairoPNG("mclu00.png", 800, 800)
# put plot command here
hist(simerups, breaks=20, probability=TRUE, 
     main="Simulated Eruption Durations", xlab="Eruption Duration (min)", col="lightgreen")

# Overlay kernel density estimate
lines(density(simerups), col="red", lwd=2)
dev.off()

# To  compare there two you would
# summary(faithful$eruptions)
# summary(simerups)
# 
# sd(faithful$eruptions)
# sd(simerups)

# by those numbers they're pretty close, though the graphs are not so convincing.

