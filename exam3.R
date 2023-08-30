#!/usr/bin/env Rscript
# this script does what?
# from chatgpt
# "r code a model of a class of students of different abilities sitting an exam with 100 questions"
library(ggplot2)
library(Cairo)

# Set seed for reproducibility
set.seed(123)

# Number of students in the class
nstu <- 50

# Number of exam questions
nques <- 100

# True underlying mean and standard deviation for student abilities
true_mean_ability <- 70
true_sd_ability <- 10

# Generate normally distributed abilities for each student
student_abilities0 <- rnorm(nstu, mean = true_mean_ability, sd = true_sd_ability)

CairoPNG("exam30.png", 800, 800)
gg <- ggplot(data = data.frame(x = student_abilities0), aes(x = x)) +
  geom_density(fill = "skyblue", color = "black") +
  labs(title = "Density Plot", x = "Data Points") +
  theme_minimal()
show(gg)
dev.off()

# Clip abilities to ensure they are within a reasonable range (optional)
student_abilities <- pmax(pmin(student_abilities0, 100), 0)

# Probability of answering each question correctly based on student ability
# You can adjust this function to reflect the relationship between ability and probability of correct response
probability_correct <- function(ability) {
  # Example: Probability increases linearly with ability
  prob_min <- 0.3
  prob_max <- 0.9
  slope <- (prob_max - prob_min) / (100 - 0)
  intercept <- prob_min
  return(pmin(pmax(intercept + slope * ability, prob_min), prob_max))
}

# Generate exam data for each student based on their abilities
pcorr <- probability_correct(student_abilities)
exam_data <- matrix(rbinom(nstu*nques, size = 1, prob = pcorr), nrow = nstu, byrow=T)

# Convert exam data to a data frame for easier analysis
exam_df <- as.data.frame(exam_data)

# Add a column for student abilities (optional)
# exam_df$Ability <- student_abilities
