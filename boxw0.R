#!/usr/bin/env Rscript
# this script does what? box and whisper pltos, no les hagas ascos!
# library(ggplot2)
library(tidyverse)
library(Cairo)

df <- data.frame(sample1=runif(4),
                 sample2=runif(4),
                 sample3=runif(4))

rownames(df) <- c("gene1","gene2","gene3","gene4")

#         sample1    sample2   sample3
# gene1 0.7068424 0.81313273 0.1021884
# gene2 0.2212768 0.87664923 0.3599538
# gene3 0.7835704 0.08712978 0.7942733
# gene4 0.3909335 0.70202803 0.8851641
# Define your responders and non-responders

responders <- c("sample1","sample2")
nonresponders <- setdiff(colnames(df),responders)

# Filter for only gene 1 and label entries

gene1 <- df[1,] %>% gather() %>% mutate(category=ifelse(key%in%responders,"responder","nonresponder"))

# Cairo image template
CairoPNG("boxw0.png", 800, 800)
qp <- qplot(x=category, y=value, data=gene1, geom=c("boxplot","jitter"), fill=category)
show(qp)
dev.off()
