#!/usr/bin/env Rscript
# this script does what?
# reor of StOv
# https://stackoverflow.com/questions/38131596/ggplot2-geom-bar-how-to-keep-order-of-data-frame
library(ggplot2)
library(Cairo)

df <- read.csv("spoca.csv", head=F)
colnames(df) <- c("Snam", "Lnam", "Val")

# lock in factor level order
df$Lnam <- factor(df$Lnam, levels = df$Lnam)

# plot
CairoPNG("bar2.png", 800, 800)
gg <- ggplot(data=df, aes(x=Lnam, y=Val)) + 
    geom_bar(stat="identity") +
    coord_flip()
show(gg)
dev.off()
