#!/usr/bin/env Rscript
# from the lintwo.R but cleaned up without the experimental code (hopefully)
# linetw2, lines between two conditions

# One of the learnigns from this, is that ggplot trying to work from the data as much as possible
# and postprocessing can be thorny, so it's best to dress up your data frame with values and strings
# that would also look good on a graph. Meaning that they will likely be a little vernbose.
library(ggplot2)
library(Cairo)
library(gridExtra)

# fix the pseudrand seq or not?
set.seed(1)

# We're going to simulate values coming from two conditions
ns <-20 
hns <- ns/2 # half the samples
c1m <- 4 # group 1 mean
c2m <- 8 # group 2 mean
c1sd <- 1 # group 1 sd
c2sd <- 1 # group 1 sd

# d <- data.frame(G1=as.numeric(c(rnorm(hns, mean=c1m, sd=c1sd), rnorm(hns, mean=c2m, sd=c2sd))),
#                 G2=as.numeric(c(rnorm(ns, mean=c1m, sd=c1sd))),
#                 Condition=rep(paste0("C",1:2), each=hns),
#                 Pair=rep(paste0("P",1:hns), 2),
#                 Legend=rep("Samples",ns))
d <- matrix(c(rnorm(hns, mean=c1m, sd=c1sd), rnorm(hns, mean=c2m, sd=c2sd),
               rnorm(ns, mean=c1m, sd=c1sd)), nrow=ns)
d <- rbind(d, c(mean(d[1:hns,1]), mean(d[1:hns,2])))
d <- rbind(d, c(mean(d[(hns+1):ns,1]), mean(d[(hns+1):ns,2])))
colnames(d) <- paste0("G", 1:2)
df <- as.data.frame(d)
# rm(d)
df$Condition=c(rep(paste0("C",1:2), each=hns), paste0("C", 1:2))
df$Pair=c(rep(paste0("P",1:hns), 2), rep("NP", 2))
df$Legend=c(rep("Samples",ns), rep("Mean", 2))
# Those column titles somewhat verbose, but make graph come out nice.

# Introduce Legend values which are the averages for both groups.
# df <- rbind(df, c(mean(df[1:hns,1]), mean(df[1:hns,2]) , "C1", "NP", "Mean"))
# df <- rbind(df, c(mean(df[(hns+1):ns,1]), mean(df[(hns+1):ns,2]) , "C2", "NP", "Mean"))
# d <- rbind(d, c(mean(d$G2[d$Condition=="C1"]), "C1", "NP", "Mean"), c(mean(d$G2[d$Condition=="C2"]), "C2", "NP", "Mean"))
# # that rbind() seems to convert numbers to character, so let's revise the column types:
df$Condition <- factor(df$Condition)
df$Pair <- factor(df$Pair)
df$Legend <- factor(df$Legend, levels=c("Samples", "Mean"))

CairoPNG("ltw2.png", 1000, 600)
p1 <- ggplot(df[,c(1,3:5)], aes(x=Condition, y=colnames(df[,c(2,3:5)])[1], group=Pair), colour=Condition) +
        # the Group=Pair will ensure the lines between each member of a pair are correct.
        geom_point(aes(colour=Legend, shape=Legend, size=Legend)) +
        geom_path(aes(colour=Legend)) +
        # quite different to have alternate sizes for point and line, so, forget.
        # scale_color_manual(values = c("darkseagreen", "limegreen"), labels=c("Samples", "Mean")) +
        scale_color_manual(values = c("darkseagreen", "limegreen")) +
        # labs(colour = "Legend") +
        # labs or labels within an aes will cause two legends! 
        # scale_fill_manual(values = c("Samples", "Means"), name="") +
        scale_size_manual(values = c(1, 4)) +
        scale_shape_manual(values = c(1, 16)) +
        # 1 is empty circ, 16 is filled circ.
        theme_bw() +
        ggtitle(paste0(colnames(df[,c(2,3:5)])[1], " expression: ", hns, " samples, 2 conditions")) +
        theme(plot.title=element_text(size=12, hjust=0.5))
p2 <- ggplot(df[,c(2,3:5)], aes(x=Condition, y=colnames(df[,c(2,3:5)])[1], group=Pair), colour=Condition) +
        # the Group=Pair will ensure the lines between each member of a pair are correct.
        geom_point(aes(colour=Legend, shape=Legend, size=Legend)) +
        geom_path(aes(colour=Legend, size=Legend), alpha=0.3) +
        # quite different to have alternate sizes for point and line, so, forget.
        # scale_color_manual(values = c("darkseagreen", "limegreen"), labels=c("Samples", "Mean")) +
        scale_color_manual(values = c("darkseagreen", "limegreen")) +
        # labs(colour = "Legend") +
        # labs or labels within an aes will cause two legends! 
        # scale_fill_manual(values = c("Samples", "Means"), name="") +
        scale_size_manual(values = c(0.5, 2)) +
        scale_shape_manual(values = c(1, 16)) +
        # 1 is empty circ, 16 is filled circ.
        theme_bw() +
        ggtitle(paste0(colnames(df[,c(2,3:5)])[1], " expression: ", hns, " samples, 2 conditions")) +
        theme(plot.title=element_text(size=12, hjust=0.5))
lp <- list(p1, p2)
do.call("grid.arrange", c(lp, ncol=2))
dev.off()

# What won't work
# scale_y_continuous(labels=scaleFUN)
# was an effort to control precision of y values.
# you get discrete value applied to continuous scale there.

# Note colour becomes the legend in geom_point()
# that's why I just named hat variable "Legend"
