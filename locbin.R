#!/usr/bin/env Rscript
# locating binaries using Sys.which(), part of R-base
# Sys.which("plink1")
ploc <- Sys.which("plink1")
cmd <- paste(ploc, "--help", sep=" ")

outp <- system(cmd, intern=T)

cat(c(cmd, " executed.\n"))
cat("Note that R will give each line of output its own line. In effect a vector of strongs, each a line.\n")
str(outp)
