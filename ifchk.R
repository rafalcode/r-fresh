#!/usr/bin/env Rscript
# this script does what?

# check when R executes a block.
# Doubts about NULL and F.

a <- NULL
b <- FALSE
c <- F

# if(a) {
#     cat("This ain't NULL disco\n")
# } else {
#     cat("This is what NULL should do\n")
# }
# You can't do the above. if(NULL) simply causes a runtime error, the program halts.
# you can tackle this via is.null() function.

if(b) {
    cat("This ain't FALSE disco\n")
} else {
    cat("This is what FALSE should do\n")
}

if(c) {
    cat("This ain't F disco\n")
} else {
    cat("This is what F should do\n")
}
