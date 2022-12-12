#!/usr/bin/env Rscript
# from the lintwo.R but cleaned up without the experimental code (hopefully)
# generate a combined gene expression and phenotype data frame, simulated data.
set.seed(1)

# We're going to simulate values coming from two conditions
ns <-20 
hns <- ns/2 # half the samples
ng <- 16
c1m <- 8 # group 1 mean
c1sd <- 1 # group 1 sd

d <- matrix(rnorm(ng*ns, mean=c1m, sd=c1sd), byrow=T, nrow=ns)

# perturbations: do it to just two of the genes:
g2pu <- 3 # g2p gene to perturb upwards
g2pd <- 9 # g2p gene to perturb downwards
d[(hns+1):ns, g2pu] <- d[(hns+1):ns, g2pu] +8
d[(hns+1):ns, g2pd] <- d[(hns+1):ns, g2pd] -2

m1 <- c()
m2 <- c()
for(i in 1:ng) {
    m1 <- c(m1, mean(d[1:hns,i]))
    m2 <- c(m2, mean(d[(hns+1):ns,i]))
}
d <- rbind(d, m1, m2)
# d <- rbind(d, m2)
colnames(d) <- paste0("G", 1:ng)
df <- as.data.frame(d)
rm(d)
df$Condition=c(rep(paste0("C",1:2), each=hns), paste0("C", 1:2))
df$Pair=c(rep(paste0("P",1:hns), 2), rep("NP", 2))
df$Legend=c(rep("Samples",ns), rep("Mean", 2))
write.csv(df, "genexp.csv", row.names=F)
