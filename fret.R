#!/usr/bin/env Rscript
# Should one fret?
# R has pretension ( I would say) to be a functional programming language
# and in that amibtions it hides its inconsistencies I would say.
# as
# https://stackoverflow.com/questions/11738823/explicitly-calling-return-in-a-function-or-not
# says
# The function() itself returns last evaluated value even without including return() function.
# But the function call itself must be aassigned to a variable.
# a and b in both this cases do exist outside the scope of their respective functions.
# You can includea return() statement
ifoo <- function() {
    a <- 42
    b <- 114
}

jfoo <- function(a=42) {
    b <- 114
    return(a) # returns a, not b, despite fact b is the last value.
}

# Overall, you can do with out return, as you're expected to know the last value is what gets returned.
# and it may have func programming logic.
