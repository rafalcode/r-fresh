#!/usr/bin/env Rscript
# hs0.R but not with Cairo
# so it 

# Create data
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
 
# Final plot, no Cairo
X11()
par(mfrow=c(1,2))
plot(h1, col=clr1, border=F , main="Pvalue Distri in bins of .05" , xlab="Pvalue bins", xlim=c(0,1))
text(0.7, mx1*0.8, paste0("pvals<.05 = ", pct1, "%"))
plot(h2, col=clr2, border=F , main="Pvalue Distri in bins of .01" , xlab="Pvalue bins", xlim=c(0,1))
text(0.7, mx2*0.8, paste0("pvals<.01 = ", pct2, "%"))
# ' readLines(stdin())
