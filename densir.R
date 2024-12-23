#!/usr/bin/env Rscript
# Desnity pltos are a standard in R.
library(Cairo)

# we'll use the famous iris dataset of course
# 150 obs, 4 variables:
# Sepal.Length
# Sepal.Width
# Petal.Length
# Petal.Width

id <-density(iris$Petal.Width)
CairoPNG("densir0.png", 800, 800)
plot(id)
dev.off()

id <-density(iris$Petal.Length)
CairoPNG("densir1.png", 800, 800)
plot(id)
dev.off()

id <-density(iris$Sepal.Width)
CairoPNG("densir2.png", 800, 800)
plot(id)
dev.off()

id <-density(iris$Sepal.Length)
CairoPNG("densir3.png", 800, 800)
plot(id)
dev.off()

CairoPNG("distir0.png", 800, 800)
hist(iris$Petal.Width)
dev.off()

CairoPNG("distir1.png", 800, 800)
hist(iris$Petal.Length)
dev.off()

CairoPNG("distir2.png", 800, 800)
hist(iris$Sepal.Width)
dev.off()

CairoPNG("distir3.png", 800, 800)
hist(iris$Sepal.Length)
dev.off()
