# Finally I asked claude.ai to do it for me.
library(Cairo)

# Original data and density
rx <- c(-2.1, -1.3, -.4, 1.9, 5.1, 6.2)
d <- density(rx, bw = 1.3)

# Create a function to interpolate the density
target_density <- approxfun(d$x, d$y, yleft = 0, yright = 0)

# Find the maximum of the target density for the envelope
max_density <- max(d$y)
M <- max_density * 1.1  # Add 10% margin

# Proposal distribution: uniform over the range
x_min <- min(d$x)
x_max <- max(d$x)
proposal_density <- 1 / (x_max - x_min)

# Rejection sampling function
rejection_sample <- function(n, seed = NULL)
{
  if (!is.null(seed))
      set.seed(seed)
  
  samples <- numeric(n)
  accepted <- 0
  total_attempts <- 0
    
  while (accepted < n) {
      # Sample from proposal (uniform)
      x_prop <- runif(1, x_min, x_max)
    
      # Sample uniform for acceptance test
      u <- runif(1, 0, M * proposal_density)
        
      total_attempts <- total_attempts + 1
        
      # Accept if u < target density
      if (u < target_density(x_prop)) {
            accepted <- accepted + 1
            samples[accepted] <- x_prop
      }
   }
  
  accrate <- n / total_attempts
  list(samples = samples, accrate = accrate, 
                total_attempts = total_attempts)
}

# Generate samples
set.seed(123)
result <- rejection_sample(1000)

cat("Acceptance rate:", round(result$accrate * 100, 2), "%\n")
cat("Total attempts:", result$total_attempts, "\n")

# Visualization
CairoPNG("kd4.png", 1600, 800)
par(mfrow = c(2, 1), mar = c(4, 4, 3, 2))

# Plot 1: Density with proposal envelope
plot(d, main = "Target Density and Proposal Envelope", 
          xlab = "x", ylab = "Density", lwd = 2, ylim = c(0, M * proposal_density))
polygon(c(x_min, x_min, x_max, x_max), 
                c(0, M * proposal_density, M * proposal_density, 0),
                        col = rgb(1, 0, 0, 0.2), border = "red", lty = 2)
points(rx, rep(0, length(rx)), pch = 19, col = "blue", cex = 1.2)
legend("topright", 
              legend = c("Target density", "Proposal envelope", "Original data"),
                     col = c("black", "red", "blue"), 
                     lty = c(1, 2, NA), pch = c(NA, NA, 19), lwd = c(2, 1, NA))

# Plot 2: Histogram of samples vs target density
hist(result$samples, breaks = 30, probability = TRUE, 
          main = "Rejection Sampling Results", 
               xlab = "x", ylab = "Density", col = "lightblue", border = "white")
lines(d, lwd = 2, col = "red")
legend("topright", 
              legend = c("Samples", "Target density"),
                     fill = c("lightblue", NA), border = c("black", NA),
                     lty = c(NA, 1), lwd = c(NA, 2), col = c(NA, "red"))

par(mfrow = c(1, 1))
dev.off()
