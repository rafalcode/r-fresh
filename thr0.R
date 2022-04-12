a <- data.frame(sample=paste0("A", 1:50), log10_total_counts=rnorm(50, 6))

# get x points
x <- 1:length(a$log10_total_counts)
png("thr0.png")
# create plot, drop x axis
plot(x, a$log10_total_counts, ylim=c(0,8), pch=19, xaxt="n", xlab="sample")
# add dashed line and X-axis labels with rotation
abline(h=6, lty="dashed")
axis(1,1:50, labels=a$sample, las=3)

# select points to color differently
LessThan6 <- a$log10_total_counts < 6
points(x[LessThan6], a$log10_total_counts[LessThan6], pch=19, col="red")
dev.off()
