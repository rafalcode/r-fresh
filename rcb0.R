#!/usr/bin/env Rscript
# This is an exercie in RColorBrewer
# which is widely used but I think quite lowly of it
# Note if you ask for more numebrs that a certan palette name can give
# it will throw and error (  I want more than 12 from Paired - it wont give that!
# Luckily the plettes themselves are just #hex number strings.
library(ggplot2)
library(Cairo)
library(RColorBrewer)

b <- brewer.pal(12, "Paired")
b2 <- brewer.pal(8, "Accent")
b3 <- brewer.pal(9, "Set1")
b4 <- brewer.pal(12, "Set3")
b5 <- brewer.pal(8, "Dark2")
ba <- c(b, b2[6], b2[8], b3[6], b4[1], b5[6])

# Cairo image template
CairoPNG("rcb0.png", 800, 800)
image(1:17,1,as.matrix(1:17),col=ba,xlab="Paired", ylab="",xaxt="n",yaxt="n",bty="n")
dev.off()
