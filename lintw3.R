#!/usr/bin/env Rscript
# from the lintwo.R but cleaned up without the experimental code (hopefully)
# takes up from lintw2, wants to have several genes.
library(ggplot2)
library(Cairo)
library(limma)

# fix the pseudrand seq or not?
set.seed(1)

# We're going to simulate values coming from two conditions
ns <-20 
ng <- 2
hns <- ns/2 # half the samples
c1m <- 4 # group 1 mean
c2m <- 6 # group 2 mean
c1sd <- 1 # group 1 sd
c2sd <- 1 # group 1 sd

de <- rbind(c(rnorm(hns, mean=c1m, sd=c1sd), rnorm(hns, mean=c2m, sd=c2sd)),
      rnorm(ns, mean=c1m, sd=c1sd))
colnames(de) <- paste0("S", 1:ns)
rownames(de) <- paste0("G", 1:ng)

phe <- data.frame(Condition=rep(paste0("C",1:2), each=hns),
                  Pair=rep(paste0("P",1:hns), 2),
                  Legend=rep("Samples",ns))
rownames(phe) <- paste0("S", 1:ns)
# Those column titles somewhat verbose, but make graph come out nice.

# Introduce Legend values which are the averages for both groups.
# beware, the means for each condition will be added to the df.
# for ggplot convenience.
de <- rbind(de, c(mean(de[1,phe$Condition=="C1"]), "C1", "NP", "Mean"), c(mean(de$G1[phe$Condition=="C2"]), "C2", "NP", "Mean"))
# that rbind() seems to convert numbers to character, so let's revise the column types:
phe$Condition <- factor(phe$Condition)
phe$Legend <- factor(phe$Legend, levels=c("Samples", "Mean"))

# CairoPNG("ltw3.png", 800, 800)
# ggplot(d, aes(x=Condition, y=de$G1, group=Pair), colour=Condition) +
#         # the Group=Pair will ensure the lines between each member of a pair are correct.
#         geom_point(aes(colour=Legend, shape=Legend, size=Legend)) +
#         geom_path(aes(colour=Legend)) +
#         # quite different to have alternate sizes for point and line, so, forget.
#         # scale_color_manual(values = c("darkseagreen", "limegreen"), labels=c("Samples", "Mean")) +
#         scale_color_manual(values = c("darkseagreen", "limegreen")) +
#         # labs(colour = "Legend") +
#         # labs or labels within an aes will cause two legends! 
#         # scale_fill_manual(values = c("Samples", "Means"), name="") +
#         scale_size_manual(values = c(2, 4)) +
#         scale_shape_manual(values = c(1, 16)) +
#         # 1 is empty circ, 16 is filled circ.
#         theme_bw() +
#         ggtitle(paste0("Simulated expression of one gene, and ", hns, " samples in each of 2 conditions")) +
#         theme(plot.title=element_text(size=16, face="bold", hjust=0.5))
# dev.off()

# What won't work
# scale_y_continuous(labels=scaleFUN)
# was an effort to control precision of y values.
# you get discrete value applied to continuous scale there.

# Note colour becomes the legend in geom_point()
# that's why I just named hat variable "Legend"
