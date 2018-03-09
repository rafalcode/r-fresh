#!/usr/bin/env Rscript
# locating binaries using Sys.which(), part of R-base
# Sys.which("plink1")
ploc <- Sys.which("plink1")
cmd <- paste(ploc, "--help", sep=" ")

outp <- system(cmd, intern=T)

cat(c(cmd, " executed.\n"))
cat("Note that R will give each line of output its own line. In effect a vector of strings, each a line.\n")
cat("Rather than print it all out, we call str() on the output, it's be a chars vector and will say how many members/lines there are in typical R fashion, i.e.\n")
browser()
str(outp)
