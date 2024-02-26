#!/usr/bin/env Rscript
# this script does what?
library(Cairo)
library(latex2exp)

# Cairo::CairoFonts('Noto Sans')
# CairoPNG("dev20.png", 800, 800)
# Cairo::CairoFonts('Helvetica')
pdf("dev20.pdf", paper="a4")
# png("dev20.png")
# dev.new()
# x11()
plot(1:20)
# text(x=1,y=15, expression(~b^2~ "value"))
# text(x=1,y=15, TeX(sprintf(r'($\beta= %d$)', 5)))
text(x=1,y=15, bquote(beta=="5"))

# the idea with pdf is that it will create ne
plot(1:20, ann=F) # ann=F tunns of axes titles.
text(x=1,y=15, bquote(beta=="5"))
# title now, but with no main.
title(sub="However, if you specify the argument adj inside \nyour plotting function all text will\nbe adjusted. If you only want some texts \nadjusted use it inside the title function")
dev.off()
# note above that subtitle ... it actually goes to the bootom but it can only ba a line long.
# I think for captions you need ggplto2.

# dev.new(width = 550, height = 330, unit = "px")
# plot(1:15)
# # Cairo image template
# CairoPNG("fname.png", 800, 800)
# plot(1:20)
# dev.off()
