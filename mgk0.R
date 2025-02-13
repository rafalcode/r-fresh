#!/usr/bin/env Rscript
# this script does what? Using R's magick package
# Note: as often happens with R, operations without an assignation are just lost
# unless theydon't go according to plan.
# So, if they work properly, they must be assigned to a variable and then displayed in some way to see the efftect of the operation.
library(magick)

im <- image_read("Angelsatmamre-trinity-rublev-1410.jpg")
im <- image_crop(im, "200x200+50+200") # overwrite im with the cropped version.
image_write(im, path = "ic.jpg", format = "jpg")
