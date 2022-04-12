#!/usr/bin/env Rscript
# first part is the main section. There are functions at the bottom
source("samputils.R")
library(lattice)
library(mice)
library('HardyWeinberg')
library('base')
source('dissort_mat.R')
source('hardy.R')
library(doParallel)
#cl <- makeCluster(2)
#registerDoParallel(cores = 4)
cores=detectCores()
cl <- makeCluster(cores[1]-1) #not to overload your computer
registerDoParallel(cl)



args <- commandArgs(trailingOnly = TRUE)
numargs <- length(args)
enumargs <- 3 
if(numargs != enumargs) {
    print("sample population number, the number of simuations, type 'yes' for dissortative model, or 'no' if you don't want it")
    stop("Stopping right here")
}
N<-as.numeric(args[1])#pop size
n<-as.numeric(args[2])#number of sim
options(warn=-1)
diss <- as.character(args[3])


sim<-function(N){
  Nd2 <- N/2 # N divided by 2 because that will be important
  #set.seed(123) # for debugging, same random sequence each time
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
  
  #if(diss == 'yes'){
  #  fsh <- dissort(ma)
  #}
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
  noff <- 1
  
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
  if(noff>1){
    for(i in 2:noff) {
      gt2 <- c(gt2, mapply(gengt, malegts, femalegts))
    }
  }

  ma2 <- cbind(ma2, gt2)
  attrinames2 <- c(attrinames2, 'GTY')
  colnames(ma2) <- attrinames2
  
  # OK we move into 3G, pairing will be more difficult now, let's order by sex
  ma2 <- ma2[order(ma2[,2]),]
  fmd <- which(ma2[,2] == 1) # female designation array
  # female designation array
  # actually mda[1] will be the number of males in G2 ... because ma2 is now ordered on sex
  fmd <- fmd[1]-1 # index is first female, so subtraction by 1 required if we want last male.
  num_f<-nrow(ma2)-fmd 
  # index is first female, so subtraction by 1 required if we want last male.
  # to decide the next mates based on sex (only)
  
  # we permute the majority sex component looking for the permutation array for G2
  #this is where the shuffling for mating happens, mating randomness
  if(fmd >= g2N/2) { # majority of males
    pa2_m <-sample(1:fmd, fmd)
    pa2_f <-((fmd+1):g2N)
  } else {
    pa2_f <-sample((fmd+1):g2N, g2N-fmd)
    pa2_m <-(1:fmd)
  }
  
  
  
  #---3rd gen---
  
  #since it's not certain it will be 50:50 males:females, I have deleted at random depending on which sex has more
  
  if(fmd > num_f){
    want_m<-pa2_m[1:num_f]
    evr <- c(want_m, pa2_f)
    ma3<-ma2[evr,]
  } else{
    want_f<-pa2_f[1:fmd]
    evr <- c(want_f, pa2_m)
    ma3<-ma2[evr,]
    
  }
#--- dissort --- here
   
  if(diss == 'yes'){
   ma3<- dissort(ma3)
   colnames(ma3) <- attrinames
  } else{
  Nma3 <- nrow(ma3)
  dNma3 <- Nma3/2
  
  FID <- ma3[1:Nma3,1]
  # so we can cbind() to our matrix, right?
  # Well, no. We need to give the female rows the number of the male they matched with
  # so the reverse of rshuf .. actually this is a tiny bit tricky
  ma3 <- cbind(ma3, c(rev(FID),FID))
  attrinames2 <- c(attrinames2, 'MTE')
  colnames(ma3) <- attrinames2
  }
  Nma3 <- nrow(ma3)
  dNma3 <- Nma3/2
  
  offnum <- ma3[1:dNma3,1] # FIDs assigned sequentially along
  ma4 <- matrix(offnum, ncol=1)
  attrinames2 <- 'FID'
  colnames(ma4) <- attrinames2
  
  # OK, what's next? Well , their sex. That's just random
  # I've turned it into a function, though it hardly needs it
  g2N <- noff*dNma3
  ma4 <- cbind(ma4, givesx(g2N))
  attrinames2 <- c(attrinames2, 'SEX')
  colnames(ma4) <- attrinames2
  
  malegts <-ma3[1:dNma3,3]
  femalegts <- rev(ma3[(dNma3+1):Nma3,3])#ma3[ma3[1:dNma3,4],3]
  gt2 <- mapply(gengt, malegts, femalegts)
  # then concatenate the rest
  if(noff>2){
    for(i in 2:noff) {
      gt2 <- c(gt2, mapply(gengt, malegts, femalegts))
    }
  }
  
  ma4 <- cbind(ma4, gt2)
  attrinames2 <- c(attrinames2, 'GTY')
  colnames(ma4) <- attrinames2
  
  #write.table(ma,file="gen1.txt")
  #write.table(ma2,file="gen2.txt")
  #write.table(ma4,file="gen3.txt")
  
  #return(list( gen1 <- data.frame(ma),gen2 <- data.frame(ma2),gen3 <- data.frame(ma4)))
  
 # gen1<-data.frame(ma)
#  gen2<-data.frame(ma2)
#  gen3<-data.frame(ma4)
  
  gty1<-ma[,3]
  gty2<-ma2[,3]
  gty3<-ma4[,3]
  
  return(list(gty1,gty2,gty3))

}



sort_rows<-function(gty){
  l0<-length(which(gty==0))
  l1<-length(which(gty==1))
  l2<-length(which(gty==2))
  return(c(l0,l1,l2))
}

Result <-list()

for(gtion in 1:3){
#gt<-matrix(nrow = n,ncol=3)
#type <- c("MM","MN","NN")
#colnames(gt) <- type
#res<-matrix(nrow=n,ncol=3)

  gt <- foreach(i = 1:n, .combine = rbind)%dopar%{
    sort_rows(sim(N)[[gtion]])
    
  
  }
  gt <- matrix(gt,ncol=3)
  res <- apply(gt,1,HWChisq)
  resher <- apply(gt,1,hetr)

  resy <- foreach(i=1:n, .combine = rbind)%dopar%{
    c(res[[i]]$pval,res[[i]]$chisq,resher[i])
  #res_type <- c("pval","chisq","heterozygosity")
  #colnames(resy) <- res_type
  }
  Result[[gtion]]<- cbind(gt,resy)  
}


stopCluster(cl)







Result<-list("Genration1" = Result[[1]],"Genration2" = Result[[2]] ,"Genration3" = Result[[3]] )

for(i in 1:3) {
  write.table(Result[[i]],file=paste0("Random_sim_gen",i,'_',N,'_',n,'_',diss,'.txt'))
}

save(Result,file = paste0('Ran_sim_result',N,'_',n,'_',diss,'.Rdata'))
