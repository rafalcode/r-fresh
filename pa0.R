#!/usr/bin/env Rscript
# I don't know who this is supposed to work.
# R uses double backslashes  .. a bit like awk.
# so this is a script which demonstrates how not to do it.

cw <- getwd()
cw
gsub("(.+)/.+$", "$1", cw, perl=T)
cw
