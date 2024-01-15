#!/usr/bin/env Rscript
# this script does what?
# yet another RoC attempt
# from https://changjunlee.com/blogs/posts/4_confusion_mat_and_roc

library(titanic)
library(Cairo)
# BTW nothing is ls() after this
library(dplyr)
library(caret)
library(pROC)

set.seed(42)

# Load the Titanic dataset
data("titanic_train") # instantiates titanic_train variable: 891 obs of 12 vars.
data("titanic_test") # 418 obs of 11 vars. It's missing the survived column, because it's a test dset.

# Combine above two:
titanic_data <- dplyr::bind_rows(titanic_train, titanic_test) # titanic_test obs get a Sruvived col, but vals are all NAs.
# 
# Select relevant features and remove rows with missing values
titanic_data <- titanic_data[, c("Survived", "Pclass", "Sex", "Age", "SibSp", "Parch", "Fare")]
titanic_data <- na.omit(titanic_data)
# so that last step is a little useless.
# any row with an NA is omitted ... so naturally all the test rows will be omitted! After having just appended them!
# 
# # Convert the 'Sex' variable to a factor
titanic_data$Sex <- as.factor(titanic_data$Sex)

train_index <- caret::createDataPartition(titanic_data$Survived, p = 0.8, list = FALSE)
titan_train <- titanic_data[train_index, ] # 
titan_test <- titanic_data[-train_index, ] # the Survived is included here.

# Create the logistic regression model
model <- stats::glm(Survived ~ ., data = titan_train, family = "binomial") # glm, as lm(), i sin stats packages. family special to glm()
# ~ . i.e. all vars except Survived are predictors

# Make predictions on the test dataset
predicted_probs <- stats::predict(model, titan_test, type = "response") # each row in test set gets a probability
predicted_classes <- ifelse(predicted_probs > 0.5, 1, 0)

# Create the confusion matrix
cm <- caret::confusionMatrix(table(predicted_classes, titan_test$Survived))

# Create the ROC curve: one function does all the work.
roc_obj <- pROC::roc(titan_test$Survived, predicted_probs)

# Setting levels: control = 0, case = 1
# Setting direction: controls < cases

# Plot the ROC curve
CairoPNG("roc2.png", 800, 800)
plot(roc_obj, main = "ROC Curve for the Logistic Regression Model")
abline(0, 1, lty = 2, col = "gray")  # Add a reference line for a random classifier
dev.off()
