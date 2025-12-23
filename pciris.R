# thi sis a simple PCA exercise, getting hte x compomenet from prcomp
# rathe rthan depending on ggplot's autoplot.
library(ggplot2)

# options(scipen=999)  # turn-off scientific notation like 1e+48
theme_set(theme_bw())  # pre-set the bw theme.

# subsampling iris a bit
df <- iris[c(10:19, 60:69, 110:119),]
pca <- prcomp(df[1:4])
dfpc <- as.data.frame(pca$x) # better off this way


# Scatterplot
gg <- ggplot(dfpc$x, aes(x=PC1, y=PC2)) + 
    geom_point(aes(col=df$Species)) # actually, yes, follows order.

#   geom_point(aes(col=state, size=popdensity)) + 
#   geom_smooth(method="loess", se=F) + 
#   xlim(c(0, 0.1)) + 
#   ylim(c(0, 500000)) + 
#   labs(subtitle="Area Vs Population", 
#        y="Population", 
#        x="Area", 
#        title="Scatterplot", 
#        caption = "Source: midwest")

# png("gg0.png")
plot(gg)
# dev.off()
