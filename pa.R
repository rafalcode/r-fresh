#!/usr/bin/env Rscript

cdup3 <- function() {
    cwd <- getwd()
    return cwd
}
cw <- cdup3()
cw
