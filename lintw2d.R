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

plotdeg0 <- function(i, df, ng)
{
    pertcols <- c(i, ng+1, ng+2, ng+3) # the pertinent columns
    print(colnames(df[,pertcols])[1])
    gg <- ggplot(df[,pertcols], aes_string(x="Condition", y=colnames(df[,pertcols])[1], group="Pair"), colour="Condition") +
          # the Group=Pair will ensure the lines between each member of a pair are correct.
          geom_point(aes_string(colour="Legend", shape="Legend", size="Legend")) +
          geom_path(aes_string(colour="Legend")) +
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
          ggtitle(paste0(colnames(df[,pertcols])[1], " expression: ", hns, " samples, 2 conditions")) +
          theme(plot.title=element_text(size=12, hjust=0.5))
      return(gg)
}

# We're going to simulate values coming from two conditions
ns <-20 
hns <- ns/2 # half the samples
ng <- 9
c1m <- 4 # group 1 mean
c2m <- 8 # group 2 mean
c1sd <- 1 # group 1 sd
c2sd <- 1 # group 1 sd

d <- matrix(rnorm(ng*ns, mean=c1m, sd=c1sd), byrow=T, nrow=ns)
m1 <- c()
m2 <- c()
for(i in 1:ng) {
    m1 <- c(m1, mean(d[1:hns,i]))
    m2 <- c(m2, mean(d[(hns+1):ns,i]))
}
d <- rbind(d, m1, m2)
# d <- rbind(d, m2)
colnames(d) <- paste0("G", 1:ng)
df <- as.data.frame(d)
rm(d)
df$Condition=c(rep(paste0("C",1:2), each=hns), paste0("C", 1:2))
df$Pair=c(rep(paste0("P",1:hns), 2), rep("NP", 2))
df$Legend=c(rep("Samples",ns), rep("Mean", 2))
write.csv(df, "genexp.csv", row.names=F)
csv <- read.csv("genexp.csv")
# set right types:
df$Condition <- factor(df$Condition)
df$Pair <- factor(df$Pair)
df$Legend <- factor(df$Legend, levels=c("Samples", "Mean"))

CairoPNG("ltw2d.png", 1000, 1000)
lp <- lapply(1:ng, FUN=plotdeg0, df=df, ng=ng)
# do.call("grid.arrange", c(lp, nrow=2, ncol=2))
do.call("grid.arrange", c(lp, ncol=3))
dev.off()
