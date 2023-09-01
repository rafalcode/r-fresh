#!/usr/bin/env Rscript
# this script does what?
# afiarly naff badly mispelt tute at
# https://www.kaggle.com/code/saife245/reticulate-a-journey-of-python-in-r

library(ggplot2)
library(Cairo)

# Cairo image template
# CairoPNG("fname.png", 800, 800)
# put plot command here
# dev.off()

library(knitr)
library(reticulate)

# knitr::knit_engines$set(python3 = reticulate::eng_python)


Sys.which("python") # thi sis R
os <- import("os") # import is from reticulate.
# os$listdir(".") # this is how we invoke a python module function 

# are we getting the right python?
use_virtualenv("r-reticulate")
sns <- import('seaborn')
data <- sns$load_dataset("fmri") # this is dataset within seaborn, becomes an R dataframe!

matplotlib <- import("matplotlib")
# matplotlib$use("Agg", force = TRUE)
plt <- import('matplotlib.pyplot')
pd <- import('pandas')

#using R's inbuilt AirPassengers dataset
df <- datasets::AirPassengers


#converting Time-Series object into an R Dataframe 
#Thx: https://stackoverflow.com/questions/5331901/transforming-a-time-series-into-a-data-frame-and-back
df1 <- data.frame(tapply(df, list(year = floor(time(df)), month = month.abb[cycle(df)]), c))
df1 <- df1[month.abb]


#building a heatmap using seaborn 
#please note the function r_to_py() that converts R object into a python 
sns$heatmap(r_to_py(df1), fmt="g", cmap ='viridis')

#display the plot
plt$show()

# I'd dearly like to know what image viewer comes up here
