# Note that the resulting graphic will change for each run,
# as the result is dependant on randomly selected values.
#
# uncomment to generate the same results with every run
set.seed(1)

# "Tight" x-axis (so it appears all horizontal lines span infinitely both left and right
par(xaxs="i")

# Sequence between -3.1 and 3.1 with 0.1 steps
x <- seq(-3.1, 3.1, 0.1)

# Plot a perfect bell curve, or a normal distribution probability density function
plot(x, dnorm(x), "l", ylim=c(0, 0.5), col="grey", lwd=5)

# 100 normally distributed random values with mean 0 and standard deviation 1
rx <- rnorm(100)

# Add lines that attempt to characterize the distribution probability density for "rx"
# using different bandwidth "bw" values 0.3, 0.1 and 0.5 (coarse to fine widths)
lines(density(rx, bw=0.3), col="#333399")
lines(density(rx, bw=0.1), col="#009900")
lines(density(rx, bw=0.05), col="#FF6633")

# Ticks at bottom to indicate where the random "rx" values are; add legend
rug(rx)
legend("topleft", legend=c("reference", "0.3", "0.1", "0.05"), bty="n",
       col=c("grey", "#333399", "#009900", "#FF6633"), lwd=c(5, 1, 1, 1))
