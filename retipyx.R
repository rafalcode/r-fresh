#!/usr/bin/env Rscript
# this script does what?
#ref. https://stackoverflow.com/questions/58982614/how-to-use-rs-reticulate-package-alongside-pythons-openpyxl-to-hide-rows-in-ex
#
library(reticulate)
library(magrittr)

rc1 <- letters
rc2 <- seq(1,1040,40)
rc3 <- seq(0,18199,700)

rc <- data.frame(rc1,rc2,rc3)

pyxl <- import("openpyxl")
wbk <- pyxl$Workbook()
pcta <- wbk$active
pcta$title <- "Trial"

for (i in 1:length(rc1)) {
    for (j in 1:length(rc)) {
        a <- rc[i,j]
        a %>% pcta$cell(i,j,.)
    }
}
# his first attempt to hide rows:
# pcta$row_dimensions[1:4] <- 10
# his second successful:
# pcta$row_dimensions$group(as.integer(1),as.integer(5),hidden = TRUE)
wbk$save("retipyx.xlsx")
