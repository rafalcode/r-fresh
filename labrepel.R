#!/usr/bin/env Rscript
# this script does what? Experiemtning to see if I can add more text to labels
# which pint out dots from afar.
library(ggrepel)
library(stringr)
library(Cairo)

set.seed(42)

dat2 <- subset(mtcars, wt > 3 & wt < 4)

# Create a new attribute called car (yes, there wasn't one before!)
dat2$car <- ""

# Let's just label these items.
ix_label <- c(2, 3, 14)
# dat2$car[ix_label] <- stringr::str_wrap(paste0(rownames(dat2)[ix_label]," hp=", dat2$hp[ix_label]), 20)
# dat2$car[ix_label] <- paste0(bquote(~bold(rownames(dat2)[ix_label])~),"\nhp=", dat2$hp[ix_label])
# dat2$car[ix_label] <- paste0(bquote(bold(.(rownames(dat2)[ix_label]))),"\nhp=", dat2$hp[ix_label])

dat2$car[ix_label] <- paste0(rownames(dat2)[ix_label],"\nhp=", dat2$hp[ix_label], "\ncyl=", dat2$cyl[ix_label])
# dat2$car[ix_label] <- bquote( "This " .(rownames(dat2)[ix_label]), splice=T)

dat3 <- rbind(
# create new dataframe of alot of extra unnamed cars
  data.frame(
    wt  = rnorm(n = 500, mean = 3),
    mpg = rnorm(n = 500, mean = 19),
    car = ""
  ),
  dat2[,c("wt", "mpg", "car")]
)

gg <- ggplot(dat3, aes(wt, mpg, label = car)) +
  geom_point(data = dat3[dat3$car == "",], color = "grey50") +
  # note here that lineheight WILL control interline distance
  # hjust=0 left alligns everything.
  # geom_text_repel(lineheight=.75, hjust=0, vjust=3, size=6, box.padding = 0.5,
  #                 max.overlaps = Inf, min.segment.length=4) +
  geom_label_repel(lineheight=.8, hjust=0, vjust=.5, size=6, # box.padding = 0.5,
                   label.padding=.25, fill="white",
                   max.overlaps = Inf, min.segment.length=4) +
  geom_point(data = dat3[dat3$car != "",], color = "red")
CairoPNG("labrepel.png", 800, 800)
show(gg)
dev.off()
