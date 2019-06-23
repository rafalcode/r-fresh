#!/usr/bin/env Rscript
fcn <- function(method){
    if(class(method) != "character") {
        stop(paste0("This function only takes a single character argument.\n"))
    }
    if (!(method == "regress" | method == "median" | method == "ks.test")) {
        stop(paste0("Error. Method \"", method, "\" is not known to this function.\n"))
    }
    return(nchar(method))
}
