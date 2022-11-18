#!/usr/bin/env Rscript
# from kaggle: https://www.kaggle.com/code/milesh1/receiver-operating-characteristic-roc-curve-in-r/notebook
# it seems you can get ROC's from randomForests as well.
# doesn't seem overly concerned with ROCs and deciding thresholds.
# Note this is RED vinho verde, which is not interesting at all.
library(Cairo)
library(randomForest)
library(pROC)

df <- read.csv("winequality-red.csv")

# change the quality column to a 1 (if quality > 5) or a 0
df$quality <- ifelse(df$quality>5,1,0)
# and then to factor
df$quality <- as.factor(df$quality)

# split the dataframe into train (80%) and test sets
index <- sample(1:nrow(df),size = 0.8*nrow(df)) # randomly selected.
train <- df[index,]
test <- df[-index,]

# build the random forest model and test it
rf_model <- randomForest(quality ~., data = train)
rf_prediction <- predict(rf_model, test, type = "prob")

# build the logistic regression model and test it
lr_model <- glm(quality ~., data = train, family = "binomial")
lr_prediction <- predict(lr_model, test, type = "response") # bit different from rf predict .. only one column.

# ROC curves
ROC_rf <- roc(test$quality, rf_prediction[,2])
ROC_lr <- roc(test$quality, lr_prediction)
# 
# Area Under Curve (AUC) for each ROC curve (higher -> better)
ROC_rf_auc <- auc(ROC_rf)
ROC_lr_auc <- auc(ROC_lr)

# plot ROC curves
CairoPNG("wineroc.png", 800, 800)
plot(ROC_rf, col = "green", main = "ROC For Random Forest (GREEN) vs Logistic Regression (RED)")
lines(ROC_lr, col = "red")
dev.off()
# print the performance of each model
print(paste0("Accuracy % of random forest: ", mean(test$quality == round(rf_prediction[,2], digits = 0))))
print(paste0("Accuracy % of logistic regression: ", mean(test$quality == round(lr_prediction, digits = 0))))
print(paste0("Area under curve of random forest: ", ROC_rf_auc))
print(paste0("Area under curve of logistic regression: ", ROC_lr_auc))
