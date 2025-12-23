# In relation to mixed effects models.
# CLaude offer this

# Mixed Effects Models in R
# Using the lme4 package

# Install if needed: install.packages("lme4")
library(lme4)
library(ggplot2)

# Example 1: Simple Random Intercept Model
# ==========================================
# Scenario: Testing a drug's effect on blood pressure
# Multiple measurements per patient over time

set.seed(123)
n_patients <- 30
n_measurements <- 5

# Simulate data
data <- data.frame(
  patient_id = rep(1:n_patients, each = n_measurements),
  time = rep(1:n_measurements, times = n_patients),
  treatment = rep(c("Drug", "Placebo"), each = n_measurements * n_patients / 2)
)

# Add patient-specific random intercepts
patient_effects <- rnorm(n_patients, mean = 0, sd = 10)
data$patient_effect <- patient_effects[data$patient_id]

# Generate blood pressure with treatment effect
data$bp <- 120 + # baseline
  -5 * (data$treatment == "Drug") + # drug effect
  -0.5 * data$time + # time trend
  data$patient_effect + # patient variation
  rnorm(nrow(data), 0, 5) # residual noise

# Fit mixed effects model
model1 <- lmer(bp ~ treatment + time + (1 | patient_id), data = data)

# View results
summary(model1)

# Extract fixed effects
fixef(model1)

# Extract random effects (patient-specific intercepts)
head(ranef(model1)$patient_id)

# Compare to naive model (ignoring clustering)
naive_model <- lm(bp ~ treatment + time, data = data)
summary(naive_model)


# Example 2: Random Slope and Intercept Model
# =============================================
# Different patients may respond differently to time

model2 <- lmer(bp ~ treatment + time + (1 + time | patient_id), data = data)
summary(model2)


# Example 3: Nested Data Structure
# ==================================
# Students within classrooms within schools

# Simulate educational data
n_schools <- 10
n_classes_per_school <- 3
n_students_per_class <- 20

edu_data <- expand.grid(
  student_id = 1:n_students_per_class,
  class_id = 1:n_classes_per_school,
  school_id = 1:n_schools
)

edu_data$student_id <- 1:nrow(edu_data)
edu_data$teaching_method <- rep(c("Traditional", "Interactive"),
                                 length.out = nrow(edu_data))

# Generate test scores with nested random effects
school_effects <- rnorm(n_schools, 0, 8)
class_effects <- rnorm(n_schools * n_classes_per_school, 0, 5)

edu_data$test_score <- 70 + # baseline
  10 * (edu_data$teaching_method == "Interactive") + # method effect
  school_effects[edu_data$school_id] + # school variation
  class_effects[edu_data$class_id + (edu_data$school_id - 1) * n_classes_per_school] +
  rnorm(nrow(edu_data), 0, 10) # student variation

# Fit nested model
model3 <- lmer(test_score ~ teaching_method + (1 | school_id/class_id),
               data = edu_data)
summary(model3)


# Visualization
# =============

# Plot individual patient trajectories
ggplot(data, aes(x = time, y = bp, group = patient_id, color = treatment)) +
  geom_line(alpha = 0.5) +
  geom_smooth(aes(group = treatment), method = "lm", se = TRUE, size = 1.5) +
  labs(title = "Blood Pressure Over Time by Treatment",
       x = "Time", y = "Blood Pressure") +
  theme_minimal()

# Model diagnostics
plot(model1)
qqnorm(resid(model1))
qqline(resid(model1))


# Model Comparison
# ================
# Compare models with different random effects structures

model_simple <- lmer(bp ~ treatment + time + (1 | patient_id), data = data)
model_complex <- lmer(bp ~ treatment + time + (1 + time | patient_id), data = data)

# Use AIC/BIC for comparison
AIC(model_simple, model_complex)
BIC(model_simple, model_complex)

# Likelihood ratio test
anova(model_simple, model_complex)
