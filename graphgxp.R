#!/usr/bin/env Rscript
# from the lintwo.R but cleaned up without the experimental code (hopefully)
# linetw2, lines between two conditions

# One of the learnigns from this, is that ggplot trying to work from the data as much as possible
# and postprocessing can be thorny, so it's best to dress up your data frame with values and strings
# that would also look good on a graph. Meaning that they will likely be a little vernbose.
library(ggplot2)
library(Cairo)
library(gridExtra)
library(limma)

plotdeg0 <- function(i, df, ng, ymax, ns)
{
    pertcols <- c(i, ng+1, ng+2, ng+3) # the pertinent columns
    gg <- ggplot(df[,pertcols], aes_string(x="Condition", y=colnames(df[,pertcols])[1], group="Pair"), colour="Condition") +
          # the Group=Pair will ensure the lines between each member of a pair are correct.
          geom_point(aes_string(colour="Legend", shape="Legend", size="Legend")) +
          geom_path(aes_string(colour="Legend", size="Legend"), alpha=0.75) +
          # quite different to have alternate sizes for point and line, so, forget.
          # scale_color_manual(values = c("darkseagreen", "limegreen"), labels=c("Samples", "Mean")) +
          scale_color_manual(values = c("darkseagreen", "darkgreen")) +
          # labs(colour = "Legend") +
          # labs or labels within an aes will cause two legends! 
          # scale_fill_manual(values = c("Samples", "Means"), name="") +
          scale_size_manual(values = c(0.25, 1)) +
          scale_shape_manual(values = c(1, 16)) +
          # 1 is empty circ, 16 is filled circ.
          ylim(0,ymax) +
          theme_bw() +
          ggtitle(paste0(colnames(df[,pertcols])[1], " expression: ", ns, " samples in 2 conditions")) +
          theme(plot.title=element_text(size=12, hjust=0.5))
      return(gg)
}

# We're going to simulate values coming from two conditions
df <- read.csv("genexp.csv")

ns <-dim(df)[1]-2
ng <-dim(df)[2]-3 # becaus we've icbind()'d pheno columns 
ncols <- sqrt(ng)
ymax <- max(df[1:ns,1:ng])
hns <- ns/2 # half the samples

# set right types:
df$Condition <- factor(df$Condition)
df$Pair <- factor(df$Pair)
df$Legend <- factor(df$Legend, levels=c("Samples", "Mean"))

CairoPNG("absolutediffs.png", 1000, 1000)
lp <- lapply(1:ng, FUN=plotdeg0, df=df, ng=ng, ymax=ymax, ns=ns)
# do.call("grid.arrange", c(lp, nrow=2, ncol=2))
do.call("grid.arrange", c(lp, ncol=ncols))
dev.off()

# Now the graphing is done, we're not interested in two rows at the end with means and the final column)
phe <- df[1:ns,(ng+1):(ncol(df)-1)]
rownames(phe) <- paste0("S", 1:ns)
data <- t(log2(df[1:ns,1:ng]))
colnames(data) <- paste0("S", 1:ns)

design <- model.matrix(~Condition, data=phe)
fit <- lmFit(data, design)
fit <- eBayes(fit)
tt<- topTable(fit,coef=2)
write.csv(tt, "quantanalysis.csv")
