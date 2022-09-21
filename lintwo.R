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
dv <- data.frame(Grp=c("G1", "G2"), Mean=c(mean(d$Val[d$Grp=="G1"]), mean(d$Val[d$Grp=="G2"])))

# CairoPNG("ltw0.png", 800, 800)
# ggplot(d, aes(x=Grp, y=Val, group=Pair), colour=Grp) +
#         geom_point() +
#         geom_path(color="green") +
#         geom_point(d, aes(x=Grp), stat="mean", colour="red")
# dev.off()

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

CairoPNG("ltw2.png", 800, 800)
ggplot(d, aes(x=Grp, y=Val)) +
       geom_point() +
       stat_summary(geom="point", fun=mean, col="purple", size=5)
#        geom_path(color="green") +
       # geom_point(data=dv, mapping=aes(x=Grp, y=Mean))
#         geom_point(d, aes(x=Grp), stat="mean", colour="red")
dev.off()
# the above is troublesome, it came from 
# https://stackoverflow.com/questions/52217852/how-to-add-means-to-a-ggplot-geom-point-plot
# .e. data=dv, introducing adding a different data frame onto the same plot? Oops, recipe for disast
