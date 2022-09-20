data(esGolub)
gol <- esGolub[1:200,]
# remove the uneeded variable 'Sample' from the phenotypic data
esGolub$Sample <- NULL

# default NMF algorithm
res <- nmf(esGolub, 3) # 3 is number of columns/rank for W, the feature/basis matrix

# is done using the default algorithm and is seeded by the default seeding methods. These defaults are set in the package specific options ’default.algorithm’ and ’default.seed’ respectively.

# Handling the result
# inspect res, brunet is the default algo.
W <- fit(res)
## <Object of class:NMFstd>
## features: 200
## basis/rank: 3
## samples: 38
# The estimated target matrix can be retrieved via the generic method fitted, which returns a (generally big) matrix:
# I think what he means here is a recontruction of the  original matrix as NMF s only an approximation
# actually W %*% H should give you the same thing.

# V.hat <- fitted(res)
# you can inspect that.

# You can str() and summary(), and more detail if you specify target
summary(res, target=esGolub)

# If there is some prior knowledge of classes present in the data, some other measures about
# the unsupervised clustering’s performance are computed (purity, entropy, . . . ). Here we use the
# phenotypic variable Cell found in the Golub dataset, that gives the samples’ cell-types (it is a
# factor with levels: T-cell, B-cell or NA):
summary(res, class=esGolub$Cell)

# The basis matrix (i.e. matrix W or the metagenes) and the mixture coefficient matrix (i.e
# matrix H or the metagene expression profiles) are retrieved using methods basis and coef respectively:

# get matrix W
w <- basis(res)
dim(w)
# get matrix H
h <- coef(res)
dim(h)

# If one wants to keep only part of the factorization, one can directly subset on the NMF object
# on features and samples (separately or simultaneously). The result is a NMF object composed of
# the selected rows and/or columns:

# keep only the first 10 features
res.subset <- res[1:10,]
# class(res.subset)
# dim(res.subset)
# # keep only the first 10 samples
# dim(res[,1:10])

# subset both features and samples:
dim(res[1:20,1:10])

# Extracting metagene-specific features

# In general NMF matrix factors are sparse, so that the metagenes can usually be characterized by
# a relatively small set of genes. Those are determined based on their relative contribution to each metagene.
# Kim and Park (KimH2007) defined a procedure to extract the relevant genes for each metagene, based on a gene scoring schema.
# The NMF package implements this procedure in methods featureScore and extractFeature:

# only compute the scores
s <- featureScore(res)
# summary(s)
# compute the scores and characterize each metagene
s <- extractFeatures(res)
# str(s)

## List of 3
## $ : int [1:8] 39 74 2 91 167 190 103 174
## $ : int [1:13] 94 1 112 42 8 64 96 182 59 41 ...
## $ : int [1:5] 43 120 128 130 129
## - attr(*, "method")= chr "kim"

# Specifying the algorithm

# Built-in algorithms

# The NMF package package provides a number of built-in algorithms, that are listed or retrieved by
# function nmfAlgorithm. Each algorithm is identified by a unique name. The following algorithms
# are currently implemented (cf. Table 1 for more details):
# nmfAlgorithm()
# [1] "brunet"
# [7] "ls-nmf"
# 
# "KL"
# "pe-nmf"
# 
# "lee"
# "siNMF"
# 
# "Frobenius" "offset"
# "snmf/r"
# "snmf/l"
# 
# "nsNMF"

# The algorithm used to compute the NMF is specified in the third argument (method). For
# example, to use the NMF algorithm from Lee and Seung (Lee2001) based on the Frobenius
# euclidean norm, one make the following call:

# using Lee and Seung's algorithm
res <- nmf(esGolub, 3, 'lee')
# algorithm(res)

# To use the Nonsmooth NMF algorithm from (Pascual-Montano2006):

# using the Nonsmooth NMF algorithm with parameter theta=0.7
resns <- nmf(esGolub, 3, 'ns', theta=0.7)
# algorithm(resns)
# fit(resns)

## <Object of class:NMFns>
## features: 200
## basis/rank: 3
## samples: 38
## theta: 0.7
# Or to use the PE-NMF algorithm from (Zhang2008):
if(requireNamespace("Biobase", quietly=TRUE)){
# using the PE-NMF algorithm with parameters alpha=0.01, beta=1
reszh <- nmf(esGolub, 3, 'pe', alpha=0.01, beta=1)

# Custom algorithms also
# stop here. the vignette has more.
