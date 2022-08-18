#!/usr/bin/env Rscript
# checks out zero inflation, especially in zero inflated negbin parameters
# quite a good link here.
# In summary the term zero inflation means that "your dataset's excessive zeros will not be inflated"
# so is a method to counter effect of many zeros in the data.
# The idea is to separately run logit regression on the zero and normal negbin regression on the rest.
# ref. https://stats.oarc.ucla.edu/r/dae/zinb/
library(ggplot2)
library(Cairo)

# Wow, that link is so complicated.
# Have no time ot delve into it just for http://jtleek.com/svaseq/simulateData.html

# the context 
# Example 2. The state wildlife biologists want to model how many fish are being caught by fishermen at a state park. Visitors are asked how long they stayed, how many people were in the group, were there children in the group and how many fish were caught. Some visitors do not fish, but there is no data on whether a person fished or not. Some visitors who did fish did not catch any fish so there are excess zeros in the data because of the people that did not fish.
zinb <- read.csv("uclafish.csv")
zinb <- within(zinb, {
  nofish <- factor(nofish)
  livebait <- factor(livebait)
  camper <- factor(camper)
})
# don't like the name count for fish count:
zinb$fcount <- zinb$count
zinb$count <- NULL

# So the idea is not so much to work out how many fish are being caught from data that only indirectly covers that.
# A camper in the US is in fact a campervan. so if camper=F, probable car motorcycle.
# it may have something to do with campers implying a fishing expedition

# If you look at summar(zinb), the levels of each factor are counted out and the "counts" refeerring to #fishcaught show very many zeros (>50%) in fact ... wnd a large minority doing all the fishing. (quite close to other behaviour of shared resources)

# histplot
CairoPNG("zinb0.png", 800, 800)
ggplot(zinb, aes(fcount, fill = camper)) +
        geom_histogram() +
        scale_x_log10() +
        facet_grid(camper ~ ., margins=TRUE, scales="free_y")
dev.off()
# in that diagra you can see that the number of counts for fcount=0 is very high.
