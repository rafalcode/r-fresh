# from prcomp example
## signs are random
require(graphics)
require(Cairo)

## the variances of the variables in the
## USArrests data vary by orders of magnitude, so scaling is appropriate
p0 <- prcomp(USArrests, scale. = TRUE)
p1 <- prcomp(~ Murder + Assault + Rape, data = USArrests, scale. = TRUE)

CairoPNG("prc00.png", 800, 800)
plot(prcomp(USArrests))
dev.off()

# summary(prcomp(USArrests, scale. = TRUE))

CairoPNG("prc01.png", 800, 800)
biplot(prcomp(USArrests, scale. = TRUE))
dev.off()
