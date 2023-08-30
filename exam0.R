#!/usr/bin/env Rscript
# this script does what?
# imsulaitn gexam marks.


library(ggplot2)
library(Cairo)


# thi swas chatgpt's attempt ... it chose the easiest

# Set seed for reproducibility
set.seed(123)

# Number of data points
n <- 200

# True underlying mean and standard deviation for marks
true_mean <- 70
true_sd <- 15

# Generate normally distributed data for student marks
student_marks <- rnorm(n, mean = true_mean, sd = true_sd)

# Ensure that the marks are between 0 and 100
student_marks <- pmax(pmin(student_marks, 100), 0)

# Display the first few data points
print(head(student_marks))

# Cairo image template
CairoPNG("fname.png", 800, 800)
# put plot command here
dev.off()
