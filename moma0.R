# simple model.matrix from Matt Richie


Id = paste0("MOUSE", 1:6)
mouse <- id
age = c(1,2,3,4,5,6)
expression = age/2 + 2 + rnorm(n=6, sd=0.1)
data <- data.frame(expression, mouse, age)
# set precision of the lofatin points.
options(digits=3)

cat("Output for Figure", ploti)
mm <- model.matrix(~age) # is a structure
fit <- lm(expression~age)
