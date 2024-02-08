#!/usr/bin/env Rscript
# the horror  I had writtn in chsa1.R, but here I got somethign decent
# so there must must an error in chsa1.R
# here I use the 8 colors which are native R names
# then they're label by the letters which correspond to
# the categorical labels in the data vector
library(ComplexHeatmap)
library(Cairo)

# from circos_color.txt
# a:chartreuse4;b:coral4;c:cyan3;d:darkgoldenrod;e:brown2;f:darkmagenta;g:deeppink2;h:orangered
circol <- c("a"="chartreuse4", "b"="coral4", "c"="cyan3", "d"="darkgoldenrod", "e"="brown2", "f"="darkmagenta", "g"="deeppink2", "h"="orangered")
ha = ComplexHeatmap::HeatmapAnnotation(foo = letters[1:8],
                     col=list(foo=circol))


# So I've worked out the default height for this, 15px for Cairo
# the struct seems to metnion 5mm, 3ox per mm? That;s sounds very reasonable.
# because 96dpi is widespread and that's per insch .. it comes to .2646 mm per px
# so just less than 4 px fit in a mm. Rough, but reasonable.

# So .. what about the width, well thi sis a bit mysterious .. it's handled automatically!

CairoPNG("chsa2.png", 300, 16) # putting it one higher so you can see white strip
draw(ha)
dev.off()
