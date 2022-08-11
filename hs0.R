#!/usr/bin/env Rscript
# this script does what?
library(Cairo)

# Create data
v <- rbeta(2000, 1 , 4)
 
# Calculate histogram, but do not draw it
# don't use a single integer with breaks (numbins)
# because R will prettify behind your back.
# best to use seq()
h1 <- hist(v, seq(0,1,.05), plot=F)
h2 <- hist(v, seq(0,1,.01), plot=F)
 
# Color vector
clr1 <- ifelse(h1$breaks < .05, "Khaki", "Grey80")
clr2 <- ifelse(h2$breaks < .01, "Gold", ifelse (h2$breaks <.05, "Khaki", "Grey80"))
 
# Final plot
# Cairo image template
CairoPNG("hs0.png", 1600, 800)
par(mfrow=c(1,2))
plot(h1, col=clr1, border=F , main="Pvalue Distri in bins of .05" , xlab="Pvalue bins", xlim=c(0,1))
plot(h2, col=clr2, border=F , main="Pvalue Distri in bins of .01" , xlab="Pvalue bins", xlim=c(0,1))
dev.off()
