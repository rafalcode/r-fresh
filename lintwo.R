#!/usr/bin/env Rscript
# lintwo , lines between two conditions
# The trouble I had with this was adding the mean as a special line
# this is the zero version, probably with alot of experimentatal code.
library(ggplot2)
library(Cairo)

# fix the pseudrand seq or not?
set.seed(1)

# I don't know what this function is for
my_smry <- function(y) { 
      ysd <- sd(y, na.rm = TRUE)
  ymd <- median(y, na.rm = TRUE)
    data.frame(y = ymd, ymin = ymd - ysd, ymax = ymd + ysd)
}

# We're going to simulate values coming from two conditions
ns <-6 
hns <- ns/2 # half the samples
g1m <- 4 # group 1 mean
g2m <- 6 # group 2 mean
g1sd <- 0.5 # group 1 sd
g2sd <- 2 # group 1 sd

d <- data.frame(row.names=paste0("S", 1:ns),
                Val=c(rnorm(hns, mean=g1m, sd=g1sd), rnorm(hns, mean=g2m, sd=g2sd)),
                Grp=rep(paste0("G",1:2), each=hns),
                Pair=rep(paste0("P",1:hns), 2),
                Meta=rep(0,ns))
d <- rbind(d, c(mean(d$Val[d$Grp=="G1"]), "G1", "NP", 1), c(mean(d$Val[d$Grp=="G2"]), "G2", "NP", 1))
# mean(d$Val[d$Grp=="G2"])

# a different dataframe ... nah, abandoned.
# dv <- data.frame(Grp=c("G1", "G2"), Mean=c(mean(d$Val[d$Grp=="G1"]), mean(d$Val[d$Grp=="G2"])))

CairoPNG("ltw0.png", 800, 800)
ggplot(d, aes(x=Grp, y=Val, group=Pair), colour=Grp) +
        geom_point(aes(colour=Meta)) +
        geom_path(aes(colour=Meta))
dev.off()

# CairoPNG("ltw1.png", 800, 800)
# ggplot(d, aes(x=Grp, y=Val, group=Pair), colour=Grp) +
# ggplot(d, aes(Grp, Val)) +
#      geom_point() +
      # geom_point(d, aes(x=Val), stat="mean")
#     geom_point(stat= "summary", fun.y=mean, shape=16, size=4, color="red")
    # stat_summary(geom="pointrange", fun.data=mean_se)
#     stat_summary(aes(x=Grp, y=Val))
    # stat_summary()
     # geom_path(color="green")
#    geom_pointrange(stat = "summary", fun.data = "my_smry", color = "red")
# dev.off()

# CairoPNG("ltw2.png", 800, 800)
# ggplot(d, aes(x=Grp, y=Val)) +
#        geom_point() +
#        stat_summary(geom="point", fun=mean, col="purple", size=5)
#        geom_path(color="green") +
       # geom_point(data=dv, mapping=aes(x=Grp, y=Mean))
#         geom_point(d, aes(x=Grp), stat="mean", colour="red")
# dev.off()
# the above is troublesome, it came from 
# https://stackoverflow.com/questions/52217852/how-to-add-means-to-a-ggplot-geom-point-plot
# .e. data=dv, introducing adding a different data frame onto the same plot? Oops, recipe for disast
