#!/usr/bin/env Rscript
# this script does what?
library(ggplot2)
library(Cairo)

CairoPNG("fname.png", 800, 800)
m <- ggplot(faithful, aes(x = eruptions, y = waiting)) +
 geom_point() +
  xlim(0.5, 6) +
   ylim(40, 110)
# contour lines
m + geom_density_2d()
dev.off()

# \donttest{
# contour bands
# m + geom_density_2d_filled(alpha = 0.5)


# contour bands and contour lines
# m + geom_density_2d_filled(alpha = 0.5) +
  geom_density_2d(size = 0.25, colour = "black")

set.seed(4393)
dsmall <- diamonds[sample(nrow(diamonds), 1000), ]
d <- ggplot(dsmall, aes(x, y))
# If you map an aesthetic to a categorical variable, you will get a
# set of contours for each value of that variable
d + geom_density_2d(aes(colour = cut))

# If you draw filled contours across multiple facets, the same bins are
# used across all facets
d + geom_density_2d_filled() + facet_wrap(vars(cut))

# If you want to make sure the peak intensity is the same in each facet,
# use `contour_var = "ndensity"`.
d + geom_density_2d_filled(contour_var = "ndensity") + facet_wrap(vars(cut))

# If you want to scale intensity by the number of observations in each group,
# use `contour_var = "count"`.
d + geom_density_2d_filled(contour_var = "count") + facet_wrap(vars(cut))


# If we turn contouring off, we can use other geoms, such as tiles:
d + stat_density_2d(
  geom = "raster",
  aes(fill = after_stat(density)),
  contour = FALSE
) + scale_fill_viridis_c()

# Or points:
# not nice, but may be interesting from a retro point of view.
d + stat_density_2d(geom = "point", aes(size = after_stat(density)), n = 20, contour = FALSE)
