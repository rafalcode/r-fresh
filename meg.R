#!/usr/bin/env Rscript
# this script does what?
library(MEGENA)
library(Cairo)

# test simplest case of planar network (a 3-clique).
# the nodes would appear to be 1, 2 and 3
# a  is a vector of start nodes
# b b is a vector of finsih nodes,
# so a and b have to be the same size and describe the edges.
# w is is a weight for those edges.
a <- c(1,1,2)
b <- c(2,3,3)
w <- runif(3,0,1)

el <- cbind(a,b,w)
el <- as.data.frame(el[order(el[,3],decreasing = TRUE),])
cpfn <- calculate.PFN(edgelist = el,max.skipEdges = Inf,doPar = FALSE,num.cores = NULL)

# Note cpfn is extremely similar to el itself, except a is called row and b is called col
# quite mysterious!

cat("and now the vignette ...\n")

# input parameters
n.cores <- 24 # number of cores/threads to call for PCP
doPar <-T # do we want to parallelize?
method = "pearson" # method for correlation. either pearson or spearman. 
FDR.cutoff = 0.05 # FDR threshold to define significant correlations upon shuffling samples. 
module.pval = 0.05 # module significance p-value. Recommended is 0.05. 
hub.pval = 0.05 # connectivity significance p-value based random tetrahedral networks
cor.perm = 10 # number of permutations for calculating FDRs for all correlation pairs. 
hub.perm = 100 # number of permutations for calculating connectivity significance p-value. 

# annotation to be done on the downstream
annot.table=NULL
id.col = 1
symbol.col= 2

data(Sample_Expression) # load toy example, but the variable is called datExpr! 30 rows, 844 cols (so many sample because TCGA)

ijw <- calculate.correlation(datExpr,doPerm = cor.perm,output.corTable = FALSE, output.permFDR = FALSE)

# Raisond for row col is probab because it's a question of correlating genes with each other.
# and as WGCNA says, all the genes are connected to each other, what matters is the eweight, which is the correlation.

#### register multiple cores if needed: note that set.parallel.backend() is deprecated. 
run.par = doPar & (getDoParWorkers() == 1) 
if(run.par) {
  cl <- parallel::makeCluster(n.cores)
  registerDoParallel(cl)
  # check how many workers are there
  cat(paste0("number of cores to use:",getDoParWorkers(),"\n"))
}

el <- calculate.PFN(ijw[,1:3],doPar = doPar,num.cores = n.cores,keep.track = FALSE)
