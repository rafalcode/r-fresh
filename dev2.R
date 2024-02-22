#!/usr/bin/env Rscript
# this script does what?
library(Cairo)
library(latex2exp)

# CairoPNG("dev20.png", 800, 800)
# pdf("dev20.pdf")
png("dev20.png")
# dev.new()
# x11()
plot(1:20)
# text(x=1,y=15, expression(~b^2~ "value"))
# text(x=1,y=15, TeX(sprintf(r'($\beta= %d$)', 5)))
text(x=1,y=15, bquote(beta=="5"))
dev.off()

# dev.new(width = 550, height = 330, unit = "px")
# plot(1:15)
# # Cairo image template
# CairoPNG("fname.png", 800, 800)
# plot(1:20)
# dev.off()
