library(Cairo)
# via Gemini I was sick 19122025
# but the issue o density plots kept me entertained.
# I really have understood them properly.

# Gemini was quite frustrating with its very loose hold on terminology
# it kept saying counts/frquency when they are not equivalent at all.
# relative frequency is merely diviing each bin in the histogram with the total numebr of dta point.
# density then is mere dividing by bin width
# as the bins are most likely all the same width, it bearly changes anything - a normalisation step you could say.



# 1. Our 5 variates - this is too little 
# alhtough the idea is to "see the artefacts"
pts <- c(12, 15, 16, 18, 25)
n <- length(pts)
h <- 1.5  # Manual bandwidth choice

# 2. Create a grid of x-values to evaluate the density
x_grid <- seq(5, 35, length.out = 200)

# 3. Calculate individual 'bumps' (Kernels) for each point
# Each kernel is a normal distribution centered at the point
kernels <- sapply(pts, function(xi) dnorm(x_grid, mean = xi, sd = h))

# 4. Average the kernels to get the final density
# (Divide by n so the total area = 1)
final_kde <- rowSums(kernels) / n

# Setup plot
CairoPNG("kde0.png", 800, 800)
plot(x_grid, final_kde, type = "l", lwd = 4, col = "black",
     # ylim = c(0, max(kernels)/n * 1.5),
     ylim = c(0, max(kernels)),
     main = "KDE of 5 Variates", xlab = "Value", ylab = "Density")

# Plot each individual kernel (scaled by 1/n)
# colors <- rainbow(n)
# these are mine
colors <-c("#A6CEE3", "#1F78B4", "#B2DF8A", "#33A02C", "#FB9A99", "#E31A1C", "#FDBF6F", "#FF7F00", "#CAB2D6", "#6A3D9A", "#FFFF99", "#B15928", "#F0027F", "#666666", "#FFFF33", "#8DD3C7", "#E6AB02")
for(i in 1:n) {
  lines(x_grid, kernels[,i] / n, col = colors[i], lty = 2)
}

# Add the raw data points as a rug
rug(pts, lwd = 2, col = "red")
dev.off()
