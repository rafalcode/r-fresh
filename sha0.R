#!/usr/bin/env Rscript
# this script does what? Using the shape functions
# 
# learnings.
# quite a few  funcstion by default add to the previous plot.
# that's why emptyplot() is often needed 
# lcol in the colour of outline
# col is the fill colour, if not specified; transparent!

# with emptyplot you set the size of the plotcanvas so to speak
# though there are fairly big margin theres
library(Cairo)
library(shape)

# Cairo image template
CairoPNG("sha0.png", 800, 800)
emptyplot(xlim=c(-3,3), frame.plot=T)
plotellipse(rx=.8, ry=.3, angle=70, lcol="dodgerblue")
plotellipse(rx=.9, ry=.2, angle=20, lcol="darkorange")

#plotellipse(rx=1, ry=.6, angle=0, from=pi, to=2*pi, arrow=T, arr.pos=seq(0.1, 0.5, by = 0.1), arr.col = rainbow(5))
# five arrowheads in that one ... a half ellipse ... don't see much use ...

# plotellipse(rx=1, ry=.6, angle=30, from=pi, to=1.2*pi, col="red")
# tiny arc .. not so interesting

# plotellipse(rx=0.1, ry=0.6, from=1.5*pi, to=pi, lcol="orange", mid=c(0.2,0.2))
# plotellipse(rx = 0.1, ry = 0.6, angle = 30, from = 1.5*pi, to = pi, lcol = "orange", mid = c(0.2,0.2))

# this next one is not from shape pkg but standard R graphics. 
# so now we must use border to set the colour, not lcol
rect(-2,-2,1,2, lwd=3, border="darkgoldenrod")
dev.off()
