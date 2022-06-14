#!/usr/bin/env Rscript
# getting barplots of not counts but continuous data
# from ggplot, in order to , say, get a weighting.
# I abndon this because I wanted a gradient effetc inside the bars
# not across them
library(ggplot2)
library(Cairo)

# note that ggplot will want a dataframe by default
# and geom_bar() is NOT the way to go, but rather geom_col()!
# that's because geom_bar() by default use counts
# but I thin you can change that twith stat=somethign else.
rv <-data.frame(vals=rnorm(12, sd=0.5))
rownames(rv) <- sprintf("%02i", as.integer(rownames(rv)))

# Cairo image template
CairoPNG("bar0.png", 800, 800)
gg <- ggplot(rv, aes(x=rownames(rv), y=vals, fill=vals)) +
#         geom_col(aes(fill=vals)) +
        geom_col() +
        scale_fill_gradient2(low ="red", high = "green", midpoint=0) +
        ggtitle("Weighting applied to each sample") +
        labs(x="sample", y="surrogate variable") +
        ylim(-2,2) +
        theme(plot.title = element_text(hjust = 0.5, face="bold", size=16))

show(gg)
dev.off()

# NOTES:
# * ggtitle always justifies to left
