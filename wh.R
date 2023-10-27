#!/usr/bin/env Rscript
# this script does what? weight hieghts correlation!
# got dataset from 
# https://gist.githubusercontent.com/nstokoe/7d4717e96c21b8ad04ec91f361b000cb/raw/bf95a2e30fceb9f2ae990eac8379fc7d844a0196/weight-height.csv
library(ggplot2)
library(Cairo)

wh <- read.csv("weight-height.csv") # 10k obs

# select an try to get a lowish correlation
wh2 <- wh[1:60,]
wh3 <- wh[30:190,] # 160 ... good low r at .76
wh4 <- wh[2600:2800,]

whc3 <- cor(wh3$Height, wh3$Weight)

# Cairo image template
# CairoPNG("fname.png", 800, 800)
# put plot command here
# dev.off()
