#!/usr/bin/env Rscript
# this script does what?
# we randomly generated p-values from the beta variate, with shape 1=1, and shape2=4
# which a heavy bottom tail (most values in first third of 0-1 space).

# rbeta() is used. Actually the beta disrib is often used for the distirbution of 
# a probability ... so it seems very suited to p-values as well.
library(Cairo)

# Create data, we shall run to different graphs of this vector of beta 
# variables. 
totpts <- 2000 # total points.
v <- rbeta(totpts, 1 , 4)
 
# Calculate histogram, but do not draw it
# don't use a single integer with breaks (numbins)
# because R will prettify behind your back.
# best to use seq()
h1 <- hist(v, seq(0,1,.05), plot=F)
mx1 <-max(h1$counts)
h2 <- hist(v, seq(0,1,.01), plot=F)
mx2 <-max(h2$counts)
pct1 <- 100*h1$counts[1]/totpts
pct2 <- 100*h2$counts[1]/totpts
 
# Color vector
clr1 <- ifelse(h1$breaks < .05, "Khaki", "Grey80")
clr2 <- ifelse(h2$breaks < .01, "Gold", ifelse (h2$breaks <.05, "Khaki", "Grey80"))
 
# Final plot
# Cairo image template
CairoPNG("hs0.png", 1600, 800)
par(mfrow=c(1,2))
plot(h1, col=clr1, border=F , main="Pvalue Distri in bins of .05" , xlab="Pvalue bins", xlim=c(0,1))
text(0.7, mx1*0.8, paste0("pvals<.05 = ", pct1, "%"))
plot(h2, col=clr2, border=F , main="Pvalue Distri in bins of .01" , xlab="Pvalue bins", xlim=c(0,1))
text(0.7, mx2*0.8, paste0("pvals<.01 = ", pct2, "%"))
dev.off()
