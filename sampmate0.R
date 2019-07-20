#!/usr/bin/env Rscript
# first part is the main section. There are functions at the bottom
source("samputils.R")

# We can make this script receive arguments, but that;s not great when debugging:
# comment out fo rhte time being and hardcode the population size instead
# args <- commandArgs(trailingOnly = TRUE)
# numargs <- length(args)
# enumargs <- 1 # expected numargs
# if(numargs != enumargs) {
#     print("This script takes 1 argument: the initial population size")
#     N <- 8
#     warning("No argument was given, there for total pop size will be 20")
# } else {
#     N <- args[1]
# }
N <- 16
Nd2 <- N/2 # N divided by 2 because that will be important
set.seed(123) # for debugging, same random sequence each time
# The argument gives us the size of the entire population.
# We first have to distinguish between males and females:
# Let's create a matrix where the individuals are row, and the column says 0
# for male and 1 for females. Let's start with the first half being males and
# the second half being females.
ma <- matrix(c(1:N, rep(0,Nd2), rep(1,Nd2)), ncol=2, byrow=FALSE)
# it's worth studying this handy set up: the c() puts everything inside a 1D vector
# then the matrix() command lays it out as a Nx2 matrix. byrow=FALSE is already the default
# I put it in for clarity, it means that the function will fill out the matrix columns first
# i.e. downwards first.

# So, the first Nd2 (i.e. in R 1:Nd2) are males and the second Nd2
# are females
# But, there's nothing really interesting about them, they're just numbers right now:
# who cares which members mate with manv members? Nobody.

# So perhaps we should introduce their genotypes to make it interesting. Use an easy format:
# 0: homozygous for ref allele, 2: homozygous for alt 1: heterozygous (actually plink's addditive format)
# Right how to assign them? Well, we're only starting out, so we'll do something unsophisticated.

# Let's use the most diverse pattern in Hardy Weinberg p^2=P(0)=0.25, 2pq=P(1)=0.5 q^2=P(2)=0.25
# So for every 4 members, we have one 0, two 1's and one 2.
# Let's do that with the rep() command
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
# Nd2+1:N are the indices of the females
fshuf <- sample((Nd2+1):N, Nd2)
# so we can cbind() to our matrix, right?
# Well, no. We need to give the female rows the number of the male they matched with
# so the reverse of rshuf .. actually this is a tiny bit tricky
fshufr <- rep(0,Nd2) # create empty vector
fshufr[fshuf-Nd2] <- 1:Nd2 # this is the operation that achieves the reverse of the permutation
ma <- cbind(ma, c(fshuf, fshufr))
# Visually, this doesn't help too much ... we have Nd2 (childless!) families now
# so why not assign a number for each pair: set males 1:Nd2, so that the female family attribute is just the index of the male they paired with
ma <- cbind(ma, c(1:Nd2, fshufr))

attrinames <- c(attrinames, 'MTE', 'FID')
colnames(ma) <- attrinames

# OK, so we've randomly paired up this generation. how will the 2nd generation pan out?
# First off, how many offspring will each pair/family have?
# let's hard code that first
noff <- 2

# Let's make another matrix for this 2nd generation, right now it only has one column, the family they belong to.
# the following will assign a family ID for each of the offspring: watch it, something apparently simple that needs to be done properly!
# offnum <- as.integer(1+((1:(noff*Nd2))-1)/noff) # FIDs bunched togther
offnum <- rep(1:Nd2, noff) # FIDs assigned sequentially along

ma2 <- matrix(offnum, ncol=1)
attrinames2 <- 'FID'
colnames(ma2) <- attrinames2

# OK, what's next? Well , their sex. That's just random
# I've turned it into a function, though it hardly needs it
g2N <- noff*Nd2
ma2 <- cbind(ma2, givesx(g2N))
attrinames2 <- c(attrinames2, 'SEX')
colnames(ma2) <- attrinames2

# so,gengt() is available ... takes two single genotypes
# we could do a forloop, but R is not good at those, instead we use mapply, it allows
# two vectors to be applied to our function.
# Again a somewhat hard to read function, we only use the first 
malegts <-ma[1:Nd2,3]
femalegts <- ma[ma[1:Nd2,4],3]
# assign first mating event
gt2 <- mapply(gengt, malegts, femalegts)
# then concatenate the rest
for(i in 2:noff) {
    gt2 <- c(gt2, mapply(gengt, malegts, femalegts))
}
ma2 <- cbind(ma2, gt2)
attrinames2 <- c(attrinames2, 'GTY')
colnames(ma2) <- attrinames2

# OK we move into 3G, pairing will be more difficult now, let's order by sex
ma2 <- ma2[order(ma2[,2]),]
fmd <- which(ma2[,2] == 1) # female designation array
# actually mda[1] will be the number of males in G2 ... because ma2 is now ordered on sex
fmd <- fmd[1]-1 # index is first female, so subtraction by 1 required if we want last male.
# to decide the next mates based on sex (only)
# we permute the majority sex component looking for the permutation array for G2
if(fmd >= g2N/2) { # majority of males
    pa2 <-sample(1:fmd, fmd)
} else {
    pa2 <-sample((fmd+1):g2N, g2N-fmd)
}
