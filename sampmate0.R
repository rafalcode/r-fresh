#!/usr/bin/env Rscript
# first part is the main section. There are functions at the bottom
args <- commandArgs(trailingOnly = TRUE)
numargs <- length(args)
enumargs <- 1 # expected numargs
if(numargs != enumargs) {
    print("This script takes 1 argument: the initial population size")
    N <- 8
    warning("No argument was given, there for total pop size will be 20")
} else {
    N <- args[1]
}

# The argument gives us the size of the entire population.
# We first have to distinguish between males and females:
# Let's create a matrix where the individuals are row, and the column says 0
# for male and 1 for females. Let's start with the first half being males and
# the second half being females.
ma <- matrix(c(1:N, rep(0,N/2), rep(1,N/2)), ncol=2, byrow=FALSE)
# it's worth studying this handy set up: the c() puts everything inside a 1D vector
# then the matrix() command lays it out as a Nx2 matrix. byrow=FALSE is already the default
# I put it in for clarity, it means that the function will fill out the matrix columns first
# i.e. downwards first.

# So, the first N/2 (i.e. in R 1:N/2) are males and the second N/2
# are females
# But, there's nothing really interesting about them, they're just numbers right now:
# who care which members mate with manv members? Nobody.
# So perhaps we should introduce their genotypes to make it interesting. Use an easy format:
# 0: homozygous for ref allele, 2: homozygous for alt 1: heterozygous
# Right how to assign them? Well, we're only starting out, so we'll do something unsophisticated
# Let's use the most diverse pattern in Hardy Weinberg p^2=P(0)=0.25, 2pq=P(1)=0.5 q^2=P(2)=0.25
# So for every 4 members, we have one 0, two 1's and one 2.
# Let's do that with the rep command
ma <- cbind(ma, rep(c(0,1,1,2), N/4))
# let's give the columns names for the matrix, just to help us
# R is sometimes awkward, we need a separate string vector (attribute names) if we want to append later:
attrinames <- c('IID', 'SEX', 'GTY')
colnames(ma) <- attrinames

# OK, enough preparation, we want the mating to go ahead.
# We're going to go for pure random mating, but of course only males with females
# For the beginning we'll say everybody finds a mate (not always true)
# Actually this is the first of randmness in the program, and we're going to lean
# on a great R function called sample() ... in fact it should be called shuffle()
# but it's capable of a few extra things .. but watch what we doing now
#
# We of course need to treat the females and males separately for this so let's
# shuffle the females and assign them to males in linear order
# N/2+1:N are the indices of the females
fshuf <- sample((N/2+1):N, N/2)
# so we can cbind() to our matrix, right?
# Well, no. We need to give the female rows the number of the male they matched with
# so the reverse of rshuf .. actually this is a tiny bit tricky
fshufr <- rep(0,4) # create empty vector
fshufr[fshuf-4] <- 1:4 # this is the operation that achieves the reverse of the permutation
ma <- cbind(ma, c(fshuf, fshufr))
# actually 
# actually this is easy, just subtract N/2 from fshuf:
attrinames <- c(attrinames, 'MTE')
colnames(ma) <- attrinames

# OK, so we've randomly paired up this generation. how will the 2nd generation pan out?
# Let's make another matrix for them, let's decide on the size first one column less than ma
# how many offspring will they have? We'll hardcode that in:
ma2c <- ncol(ma)-1
noff <- 2
ma2 <- matrix(rep(0, ma2c*noff*N/2), ncol=ma2c)
ma[ma[1:4,4],3]

# gengt(): function to generate a new gt assignment given two gt's.
gengt <- function(gt1, gt2)
{
    if((gt1 == 0) && (gt2 ==0)) { 
        return(0) # two ref alleles? always homozygous for ref allele
    } else if(((gt1 == 0) && (gt2 ==2)) || ((gt1 == 2) && (gt2 ==0))) {
        return(1) # homozyg for ref and homozyg for alt? Then always heterozyg
    } else if((gt1 == 2) && (gt2 ==2)) {
        return(2) # two alt alleles? always homozygous for alt allele
    } else if(((gt1 == 0) && (gt2 ==1)) || ((gt1 == 1) && (gt2 ==0))) {
        return(round(runif(1))) # homozyg for ref and heterozyg? 50% chance of homozyg ref, 50% chance heterozyg
    } else if(((gt1 == 2) && (gt2 ==1)) || ((gt1 == 1) && (gt2 ==2))) {
        return(round(1+runif(1))) # homozyg for alt and heterozyg? 50% chance of homozyg alt, 50% chance heterozyg
    } else if((gt1 == 1) && (gt2 ==1)) {
        # both heterozyg's? slightly more complicated: 25% chance homozyg ref, 50% heterozyg, 25 homozyg alt
        r <- runif(1)
        if(r < 0.25) {
            return(0)
        } else if(r < 0.75) {
            return(1)
        } else {
            return(2)
        }
    }
}

