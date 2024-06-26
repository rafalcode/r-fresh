#!/usr/bin/env Rscript
# this script does what? From https://www.digitalocean.com/community/tutorials/plot-roc-curve-r-programming
library(Cairo)
library(caret)
library(pROC)

# https://github.com/Safa1615/Bike-loan-Dataset
data = read.csv("bank-loan_womv.csv",header=TRUE) # without missing values

### Data SAMPLING ####
set.seed(101)
split = createDataPartition(data$default, p = 0.80, list = FALSE)
train_data = data[split,]
test_data = data[-split,]

#error metrics -- Confusion Matrix
err_metric=function(CM) {
  TN =CM[1,1]
  TP =CM[2,2]
  FP =CM[1,2]
  FN =CM[2,1]
  precision =(TP)/(TP+FP)
  recall_score =(FP)/(FP+TN)
  f1_score=2*((precision*recall_score)/(precision+recall_score))
  accuracy_model  =(TP+TN)/(TP+TN+FP+FN)
  False_positive_rate =(FP)/(FP+TN)
  False_negative_rate =(FN)/(FN+TP)
  print(paste("Precision value of the model: ",round(precision,2)))
  print(paste("Accuracy of the model: ",round(accuracy_model,2)))
  print(paste("Recall value of the model: ",round(recall_score,2)))
  print(paste("False Positive rate of the model: ",round(False_positive_rate,2)))
  print(paste("False Negative rate of the model: ",round(False_negative_rate,2)))
  print(paste("f1 score of the model: ",round(f1_score,2)))
}

# 1. Logistic regression
logit_m =glm(formula = default~. ,data =train_data ,family='binomial')
summary(logit_m)
logit_P = predict(logit_m , newdata = test_data[-13] ,type = 'response' )
logit_P <- ifelse(logit_P > 0.5,1,0) # Probability check
CM= table(test_data[,13] , logit_P)
print(CM)
err_metric(CM)

#ROC-curve using pROC library
roc_score=roc(test_data[,13], logit_P) #AUC score
CairoPNG("roc0.png", 800, 800)
plot(roc_score ,main ="ROC curve -- Logistic Regression ")
dev.off()
