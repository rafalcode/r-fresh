#!/usr/bin/env Rscript
# linetwo , lines between two conditions
library(ggplot2)
library(Cairo)

# fix the pseudrand seq or not?
set.seed(1)

my_smry <- function(y) { 
      ysd <- sd(y, na.rm = TRUE)
  ymd <- median(y, na.rm = TRUE)
    data.frame(y = ymd, ymin = ymd - ysd, ymax = ymd + ysd)
}

ns <-6 
hns <- ns/2 # half the samples
g1m <- 4 # group 1 mean
g2m <- 6 # group 2 mean
g1sd <- 0.5 # group 1 sd
g2sd <- 2 # group 1 sd

d <- data.frame(row.names=paste0("S", 1:ns),
                Val=c(rnorm(hns, mean=g1m, sd=g1sd), rnorm(hns, mean=g2m, sd=g2sd)),
                Grp=rep(paste0("G",1:2), each=hns),
                Pair=rep(paste0("P",1:hns), 2))

# CairoPNG("ltw0.png", 800, 800)
# ggplot(d, aes(x=Grp, y=Val, group=Pair), colour=Grp) +
#         geom_point() +
#         geom_path(color="green") +
#         geom_point(d, aes(x=Grp), stat="mean", colour="red")
# dev.off()

CairoPNG("ltw1.png", 800, 800)
ggplot(d, aes(x=Grp, y=Val, group=Pair), colour=Grp) +
#         geom_point(d, aes(x=Val), stat="mean")
#    geom_point(stat= "summary", fun.y=mean, shape=16, size=4, color="red")
     geom_point() +
     geom_path(color="green")
#    geom_pointrange(stat = "summary", fun.data = "my_smry", color = "red")
dev.off()

