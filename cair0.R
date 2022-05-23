library(Cairo)
# Cairo(600, 600, file="cair0.png", type="png", bg="white")
CairoPDF("plot.pdf", 6, 6, bg="transparent")
plot(rnorm(4000),rnorm(4000),col="#ff000018",pch=19,cex=2) # semi-transparent red
dev.off() # creates a file "plot.png" with the above 
