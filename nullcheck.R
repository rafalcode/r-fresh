#!/usr/bin/env Rscript
# this script does what? checking out conditionals on NULL

k <- NULL
if(!is.null(k)) {
    cat("k ACTIVE\n")
} else {
    cat("k is not active\n")
}

# the lesson is that NULL is not FALSE,
# so you have to use is.null() on them.
