# from the Gemini inquiry
# the histogram object in R is quite useful.
set.seed(42)
data <- rnorm(100, mean = 170, sd = 10)

# 2. Define our bins manually
breaks <- seq(140, 200, by = 5)

# diff() will just subtract 1st from 2nd, 2nd from 3rd etc.
# it's one way of doing it.
bin_width <- diff(breaks)[1] # in this case, 5

# 3. Get the counts for each bin
counts <- hist(data, breaks = breaks, plot = FALSE)$counts

# 4. Calculate Density using first principles
# Density = (Count / Total) / BinWidth
relative_freq <- counts / length(data)
manual_density <- relative_freq / bin_width

# 5. Compare with R's built-in density calculation
r_density <- hist(data, breaks = breaks, plot = FALSE)$density
all.equal(manual_density, r_density) # Should return TRUE

