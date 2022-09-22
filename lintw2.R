#!/usr/bin/env Rscript
# from the lintwo.R but cleaned up without the experimental code (hopefully)
# linetw2, lines between two conditions
library(ggplot2)
library(Cairo)

# fix the pseudrand seq or not?
set.seed(1)
options(digits=4) # this won't work for ggplot2
# You can create function instead
scaleFUN <- function(x) sprintf("%.3f", x)

# We're going to simulate values coming from two conditions
ns <-6 
hns <- ns/2 # half the samples
c1m <- 4 # group 1 mean
c2m <- 6 # group 2 mean
c1sd <- 0.5 # group 1 sd
c2sd <- 2 # group 1 sd

d <- data.frame(row.names=paste0("S", 1:ns),
                Val=c(rnorm(hns, mean=c1m, sd=c1sd), rnorm(hns, mean=c2m, sd=c2sd)),
                Cond=rep(paste0("C",1:2), each=hns),
                Pair=rep(paste0("P",1:hns), 2),
                Meta=rep(0,ns))

# Introduce Meta values which are the averages for both groups.
d <- rbind(d, c(mean(d$Val[d$Cond=="C1"]), "C1", "NP", 1), c(mean(d$Val[d$Cond=="C2"]), "C2", "NP", 1))

CairoPNG("ltw2.png", 800, 800)
ggplot(d, aes(x=Cond, y=Val, group=Pair), colour=Cond) +
        # the Group=Pair will ensure the lines between each member of a pair are correct.
        geom_point(aes(colour=Meta, shape=Meta, size=Meta)) +
        geom_path(aes(colour=Meta)) +
        scale_color_manual(values = c("darkseagreen", "limegreen")) +
#        scale_fill_manual(values = c("Samples", "Means"), name="") +
        scale_size_manual(values = c(2, 2)) +
        scale_shape_manual(values = c(1, 16)) +
        # 1 is empty circ, 16 is filled circ.
        ggtitle(paste0("Expression of one gene, and ", hns, " samples in each of 2 conditions")) +
        guides(fill=guide_legend(title="Type")) +
        theme(plot.title=element_text(size=16, face="bold"))
dev.off()

# What won't work
# scale_y_continuous(labels=scaleFUN)
# was an effort to control precision of y values.
# you get discrete value applied to continuous scale there.

