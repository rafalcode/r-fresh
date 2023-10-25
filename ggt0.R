#!/usr/bin/env Rscript
# this script does what? Concentrating on overlaying text
library(ggplot2)
library(Cairo)
df0 <- data.frame(x = c(1, 1, 2, 2, 1.5),
                 y = c(1, 2, 1, 2, 1.5),
                 text = c("hi from bottomleft", "hi from topleft", "hi from bottomright", "hi from topright", "hi from center"))

df <- data.frame(x = c(1, 1, 12, 12, 1.5),
                 y = c(1, 12, 1, 12, 1.5))

txt = c("hi from bottomleft", "hi from topleft", "hi from bottomright", "hi from topright", "hi from center")

# ggp <- ggplot(df0, aes(x, y))
# ggp <- ggp + geom_text(aes(label = text))
ggp <- ggplot(df, aes(x, y))


minx <- min(ggp$data$x)
maxx <- max(ggp$data$x)
miny <- min(ggp$data$y)
maxy <- max(ggp$data$y)

# if you leave the following out, ggplot will evaluate the plot as being empty
# and when later annotation are added it will resize, which is nasty.
# mostly jowever I will not anotate on empty ggplot graphs.

ggp <- ggp + geom_text(aes(label = txt))
# CairoPNG("ggt0.png", 800, 800)
# show(ggp)
# dev.off()

# OK, now try and add two annotations
xs <- c(minx + (maxx-minx)*.3, minx + (maxx-minx)*.8)
ys <- c(miny + (maxy-miny)*.7, miny + (maxy-miny)*.9)
t2 <- c("onean", "twoan")

ggp <- ggp +annotate("text", x=xs, y=ys, label=t2)

CairoPNG("ggt0.png", 800, 800)
show(ggp)
dev.off()
