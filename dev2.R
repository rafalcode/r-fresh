#!/usr/bin/env Rscript
# this script does what?
library(Cairo)
library(latex2exp)

# CairoFonts("FreeSans")
# trouble with this function is that it doesn't return anything .. how to verify?

# try ...
# CairoFontMatch(fontpattern="Helvetica",sort=FALSE,verbose=F)
# 1. family: "Nimbus Sans", style: "Regular", file: "/usr/share/fonts/opentype/urw-base35/NimbusSans-Regular.otf" 

# CairoFonts("Nimbus Sans") # nothing about return value
# CairoFonts("Noto Serif", usePUA=F) # this does have an effect! definnitely changes to Serif, not sure about Noto though.
CairoFonts(regular="Noto Serif", symbol="Standard Symbols PS")
# CairoFonts("Standard Symbols PS", usePUA=F)
# CairoFonts("Comfortaa")
CairoPNG("dev20.png", 800, 800)
# pdf("dev20.pdf")
# png("dev20.png")
# dev.new()
# x11()
plot(1:20)
# text(x=1,y=15, expression(~b^2~ "value"))
# text(x=1,y=15, cex=2, TeX(sprintf(r'($\beta= %d$)', 5)))
# text(x=1,y=15, cex=2, adj=c(0,0), bquote("The greek letter beta looks like"~beta))
text(x=1,y=15, cex=2, adj=c(0,0), "The greek letter beta looks like β")
# vim to enter those bad boys
# ctrl+k, b, *
# gives beta, rho alpha similar ... what theta it's h*
dev.off()

# dev.new(width = 550, height = 330, unit = "px")
# plot(1:15)
# # Cairo image template
# CairoPNG("fname.png", 800, 800)
# plot(1:20)
# dev.off()
