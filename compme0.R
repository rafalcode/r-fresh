#!/usr/bin/env Rscript
# this is the stat_compare_means routine from the refman.
# you will get "cannot compute exact p-value due to ties" warning
# you can suppress them with exact=F ... but that doesn;t take away the problem of course

# Important note, this seems to highlight some problem Cairo and the show() function
# no this is just my silly mistakeA p <- p + etc() is required!

library(ggplot2)
library(ggpubr)
library(Cairo)

data("ToothGrowth")
# 10 Guinea pigs for 3 types of dosage and 2 types of supplements Vitamin C(VC) or Orange Jiuce (OJ).
# len is the dependent outcome variable.

# Two independent groups
png("cme0.png", 800, 800)
p0 <- ggboxplot(ToothGrowth, x = "supp", y = "len", color = "supp", palette = "npg", add = "jitter")
show(p0)
dev.off()

# Add p-value
CairoPNG("cme1.png", 800, 800)
p <- p0 + stat_compare_means()
show(p)
dev.off()
# just no difference between those two.
# paired=T does nothing.

# Default method is Wilcoxon.
CairoPNG("cme2.png", 800, 800)
p <- p0 + stat_compare_means(method = "t.test", label="p.format")
show(p)
dev.off()
# No difference in the graphs!

# Paired samples
CairoPNG("cme3.png", 800, 800)
p <- ggpaired(ToothGrowth, x = "supp", y = "len",
             color = "supp", line.color = "gray", line.size = 0.4, palette = "npg") +
    stat_compare_means(paired = TRUE)
show(p)
dev.off()

# More than two groups
# Pairwise comparisons: Specify the comparisons you want
my_comparisons <- list( c("0.5", "1"), c("1", "2"), c("0.5", "2") )
CairoPNG("cme4.png", 800, 800)
p <- ggboxplot(ToothGrowth, x = "dose", y = "len", color = "dose", palette = "npg")+
    # Add pairwise comparisons p-value
    stat_compare_means(comparisons = my_comparisons, label.y = c(29, 35, 40))+
    stat_compare_means(label.y = 45)     # Add global Anova p-value
show(p)
dev.off()

# Multiple pairwise test against a reference group
CairoPNG("cme5.png", 800, 800)
p <- ggboxplot(ToothGrowth, x = "dose", y = "len",
              color = "dose", palette = "npg") +
    stat_compare_means(method = "anova", label.y = 40) + # Add global p-value
    stat_compare_means(aes(label = ..p.signif..), method = "t.test", ref.group = "0.5")
show(p)
dev.off()

# Multiple grouping variables
# Box plot facetted by "dose"
CairoPNG("cme6.png", 800, 800)
p <- ggboxplot(ToothGrowth, x = "supp", y = "len", color = "supp", palette = "npg",
               add = "jitter", facet.by = "dose", short.panel.labs = FALSE) +
                # Use only p.format as label. Remove method name.
                stat_compare_means( aes(label = paste0("p = ", ..p.format..)))
show(p)
dev.off()
