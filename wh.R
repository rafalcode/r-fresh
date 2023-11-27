#!/usr/bin/env Rscript
# this script does what? weight hieghts correlation!
# got dataset from 
# https://gist.githubusercontent.com/nstokoe/7d4717e96c21b8ad04ec91f361b000cb/raw/bf95a2e30fceb9f2ae990eac8379fc7d844a0196/weight-height.csv
library(ggplot2)
library(Cairo)
library(DescTools)

wh <- read.csv("weight-height.csv") # 10k obs

# select an try to get a lowish correlation
wh2 <- wh[1:60,]
wh3 <- wh[30:190,] # 160 ... good low r at .76
wh4 <- wh[2600:2800,]

rho3 <- cor(wh3$Height, wh3$Weight) # you only get rho back from this function
# using DescTools' FIsherZ ... it doesn't need n, just the rho .76 will give a
# z of  about .99 ... z goes over 1 ... .86 is 1.29.

# Cairo image template
# CairoPNG("fname.png", 800, 800)
# put plot command here
# dev.off()
