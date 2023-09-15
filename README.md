# r-fresh

Refreshing R skills.


# pretty printing
You use cat(mystring)
no new line is added.
To get printf then
you need sprintf

# bexp.R
run a bash command via R. Critically, uses a bash expansion, so not entirely trivial.

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
uses dbeta(). ./becu 5 2 is the right biased one, and 2 5 is the left-sided one.

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

# sca0.R and sca1.R
experiments with the scale() function.

# bar2.R
What geom_bar, it may reorder

# chunks.R
I've often done this in C, somehow, nto for R.
It's actually dead easy, but it's also easy to get mixed up, it's just a basic chunkifier.

# sdgene.R
Linear models often add an extra term (multiplied by a coeeficient), but how much does that affect variation?
Well only incremental is would seem. A single effect with a higher SD, say 4, would have a much bigger effect than 4
additions of sd=1 rnorm.

# R colors for C
there are 657 of these, rcol1.R was used to generate most of rcolp.h (which means rcolor proportional RGB vlaues givien in [0,1] value format). I then editted th 

# sgd series
These are the seqgendiff trials. In the beginning it held promise, but I wrote
 G for Gerard, G for garbled.
for a reason, his terminology seems weird.
* uses null for zero
* seems obsessed with thinning all the time. This may have something to do with variance reduction, but is that right? Because he's pushing it all the time
* "databased" is the idea of using existing (real) rnaseq datasets, and adding a simulated signal to it. This is the approach of seggendiff.
* thinning, actually is binomial thinning and it means subsampling from the bionmial distribution and this is actually the method of actuall applying the artificial signal. Perhaps isntead of just adding? Maybe. 

# Rmd files.
there is more than one version it appears. v2 must be run via the rmarkdown package. see rmd2h.R script. dia.Rmd was the input file for it.
