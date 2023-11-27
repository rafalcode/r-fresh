#!/usr/bin/env Rscript
# this script does what?
library(ggplot2)
library(Cairo)

# More of a general enquiry into density plots.
# I particular wonder if smoothing functions are used to make them look better
# try chatgpt
# Here's ggplot
# ggplot2 density plot
# 
# To create a density plot using ggplot2 in R, we'll first need to install the package (if you haven't done it already) and then create the plot. We'll use a sample dataset to demonstrate the density plot. In this example, let's create a density plot for a random set of data.
# 
# Load the ggplot2 package
library(ggplot2)

# Generate a random dataset with <#yourchoice> data points
howmanydps <- 10
set.seed(123)
random_data <- rnorm(howmanydps, mean = 50, sd = 10)

# Create a density plot using ggplot2
p <- ggplot(data = data.frame(x = random_data), aes(x = x)) +
      geom_density(fill = "skyblue", color = "black") +
        labs(title = "Density Plot", x = "Data Points") +
          theme_minimal()

#Chatgpt's explanation:
# In this example, we generated a random dataset of 100 data points from a normal distribution with a mean of 50 and a standard deviation of 10. Then, we used ggplot2 to create a density plot using geom_density().
# 
# The aes() function is used to specify the aesthetics, and we map the variable x (our data points) to the x-axis. We use geom_density() to add the density layer to the plot, with fill = "skyblue" specifying the color of the density area and color = "black" specifying the color of the outline.
# 
# We also set the plot title and label for the x-axis using labs(), and we use theme_minimal() to set a minimalistic theme for the plot. Of course, you can customize the plot further by adding other ggplot2 elements like axes labels, legends, etc., depending on your specific needs and preferences.


# So, is there smoothing going on?  You bet! But there is also a bit more too.
# perhaps the number bins are normalised too.
# Caire image template
CairoPNG("dens1.png", 800, 800)
show(p)
dev.off()
