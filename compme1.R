#!/usr/bin/env Rscript

# Pairwise comparisons: Specify the comparisons you want
my_comparisons <- list( c("0.5", "1"), c("1", "2"), c("0.5", "2") )
CairoPNG("cme4a.png", 800, 800)
p <- ggboxplot(ToothGrowth, x = "dose", y = "len", color = "dose", palette = "npg")+
    # Add pairwise comparisons p-value
    stat_compare_means(comparisons = my_comparisons, label.y = c(29, 35, 40)) +
    stat_compare_means(label.y = 45)     # Add global Anova p-value (Kruskal-Wallis)
    # note that the label.y are vertical positions and are mapped to the range of the y axis.
    # in this case, it's previously worked out that 45 is a good vertical level for this particular plot.
show(p)
dev.off()
