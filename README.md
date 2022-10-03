# r-fresh

Refreshing R skills.


# pretty printing
You use cat(mystring)
no new line is added.
To get printf then
you need sprintf

# rcpp
cxxfunction that's not from Rccp but from inline

# sto0.R
using stop() with a function, and also, how to escape quotations.

# fustop0.R
Further adventures in using stop(), this time showing how stop(-1) is worse than useless.

# the f2s\*
scripts are some way of debugging R when arguments are being called in 

# sta0.R and sta1.R
static variables in R. Upshot is you need a closure
of a function within a function with a <<- local assignment.

# becu.R
Printing out a beta distribution, ran into a few fundamental problems here.

# desvig0.R
Simple 2 genes interaction example from the DESeq2 vignette.

# tracol.R
grDevices:color() will get the list of 600 or R colors which actually have names.
it's handy until you want opacity, named colors won;t do that, you need a function
for it ... here it is.

# lintwo.R and lintw2.R
Diff expression of one gene, though more an exercise in ggplot.

# tangheat.R and tangheat2.R
Using the pheatmap package because it seems to be well-known for its annotation abilities.
It does have limitations, the legend positions area bit hard-coded.

# ddm
degenerate matrices.
