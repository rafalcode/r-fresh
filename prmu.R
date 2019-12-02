#!/usr/bin/env Rscript
# Script highlighting permutations.
# Generating permutation is easy
# but generating the reverse is trickier.
#
# One thing must be said: you can;t generate a reverse permutation
# without any extra information
# i.e. you have an array or vctor of integers in a certain order, right. Well you can't
# get back to its "original order" because you need extra information for that.
#
# the best way is to have have an array of indices, and to assume that the previous 
# order is the canonical monotonically increasing 1 to N
#
# You need to have you mind clear the distinction between a vectors' indices and its values at those indices.
# Also beware of the common occurence of a vector whose values are the same as the indices, this really blurs 
# the difference between values and indices, but it's not illegal either.
# 
# OK let's generate a vector where the values are quite different to the indices:
vecint <- as.integer(10+40*runif(8))
vecint
# 8 random numbers between 10 and 50.

# Now you generate a permutation very easily with sample() which will
# generate shuffled (note!) values, but actually it's better to generate
# and intermediate file and shuffle the indices instead.
#
# let's see: enter shufidx: the array of shuffled indices:
shufidx <- sample(1:length(vecint), length(vecint))
shufidx
#
# That is not our shuffled vector, but the following is
vecint[shufidx]
# Actually typicall what will happen is that you will overwrite vecint, 
# so we'd better do that:
vecint <- vecint[shufidx]
# Note that we are implicitly re-assigning to vecint in monotonically increasing indices

# now we want to recover the old vecint
# we create a zero value vector of the same size
vecintr <- rep(0, length(vecint))
# and this time we don't assign to vecintr in montonicall increasing indices
# but rather in shuffled index order. Its *values* however, are the monotonically increasing
# (newly ordered) vecint values.
vecintr[shufidx] <- vecint
vecintr
# Note that if vecint was just 1,2,3 etc. you could write
# vecintr[shufidx] <- 1:length(vecint)
# and get way with it. But it's misleading because that is the single easy case.
# Note how you have to keep the shufidx vector. That's why above I stressed that you need something more




# fshufr <- rep(0,Nd2)
# fshuf <- ma[(Nd2+1):N,1]
# fshufr[fshuf-Nd2] <- 1:Nd2
