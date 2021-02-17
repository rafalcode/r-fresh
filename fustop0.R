#!/usr/bin/env Rscript
# this s just a demo script
# to show how stop(-1) is an appalling way to exit a function when things go wrong.
#
# the message is just -1. Even a minimal string would be better, i.e.
# stop("bad string")
# because -1 is supposed to mean something and you will wate time tryign to find out.
# i.e. the system will give 1 as the return alue is there is an error. That is not -1.
# so just returning -1 is awful.
# Note that at least R will say in which function it occured. But you get that for free.

myfunc <- function(string)
{
    if( string == "doeswhat?") {
        warning("Yes you guessed string correctly")
    } else {
        # stop("Sorry, wrong string!")
        stop(-1)
    }
}

ret <- myfunc("sorry?")
cat(paste0("R script return value is ", ret, "\n"))
# myfunc("doeswhat?")
