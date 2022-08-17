#!/usr/bin/env Rscript
# checks out zero inflation, especially in zero inflated negbin parameters
# quite a good link here.
# In summary the term zero inflation means that "excessive zeros will not be inflated"
# and can prove useful if there are many zeros in the data.
# The idea is to separetely run logit regression on the zero and normal negbin regression on the 
# ref. https://stats.oarc.ucla.edu/r/dae/zinb/
library(ggplot2)
library(Cairo)


# the context 
# Example 2. The state wildlife biologists want to model how many fish are being caught by fishermen at a state park. Visitors are asked how long they stayed, how many people were in the group, were there children in the group and how many fish were caught. Some visitors do not fish, but there is no data on whether a person fished or not. Some visitors who did fish did not catch any fish so there are excess zeros in the data because of the people that did not fish.
zinb <- read.csv("uclafish.csv")
zinb <- within(zinb, {
  nofish <- factor(nofish)
  livebait <- factor(livebait)
  camper <- factor(camper)
})




# Cairo image template
# CairoPNG("fname.png", 800, 800)
# put plot command here
# dev.off()
