#!/usr/bin/env Rscript
# playing with the sample function
# although a game it also highlight the glory of the table() function
# with is hideously boringly named
# it applies to factors.
sasz <- 10 # sample size
mxpla <- 10 # max numba of playerz (actual will be lower, depending on chance)

fld <- as.factor(floor(runif(sasz) * mxpla))
tfld <- sort(table(fld), decreasing=T)
tsz <- length(tfld)
lastnz <- tsz
tfld
i <- 1
while(lastnz>1) {
    cat(paste0("Iteration #", i, ":\n"))
    fld <- sample(fld, sasz, replace=T)
    tfld <- sort(table(fld), decreasing=T)
    tsz <- length(tfld)
    retma <- match(0, tfld)
    if(is.na(retma)) {
        print(tfld[1:tsz])
    } else {
        lastnz <- retma -1 # last nonzero
        print(tfld[1:lastnz])
    }
    i <- i + 1
}

# discussion.
# I'd like to keep every iteration sort so I can print out the winning number at the end
# but I want to keep the zeros out.
# so I want to detect the first zero.

# match() is built for detecting the first element of something.
# unfort doesn't seem to have anything for the last match of something.
# so, have to sort decreasing.

# table() is slightly weird, especially when printing out 
# it always includes the name of vector it was created out of.
