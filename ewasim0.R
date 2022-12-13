library(truncnorm)
library(pwrEWAS)
library(pwrEWAS.data)

# I've checked the Colon one OK, idx 7.
tissueType = c("Adult (PBMC)",
    "Saliva", 
    "Sperm", 
    "Lymphoma",
    "Placenta",
    "Liver",
    "Colon",
    "Blood adult",
    "Blood 5 year olds",
    "Blood newborns",
    "Cord-blood (whole blood)",
    "Cord-blood (PBMC)") # tissue type that is used as reference to sample from

methPara <- pwrEWAS.data:::loadDataset(tissueType[7])
CpGonArray <- length(methPara$mu)

# J = 100000, # number of simulated CpGs
J = 10000 # very unusual variable name for the total number of CpG to simulate. Needs to less than size of methPara (270k or so).

# don't have this sorted out at all.
## sample CpGs
cpgIdx <- sample(x = seq_len(CpGonArray), size = J, replace = TRUE) # pick J random CpG's to be simulated
cpgIdxName <- paste(seq_len(J), "_", rownames(methPara)[cpgIdx], sep = "") # ensuring unique CpG name (allowing unique sampling with replacement)
stop("o")

changedCpgsIdx <- sample(x = cpgIdx, size = K[d]) # pick K random CpG's to be changed in mean meth
changedCpgsIdxName <- cpgIdxName[match(changedCpgsIdx, cpgIdx)]

## Change Mu for "changedCpgsIdx"'s CpG's
# drawing delta from truncated normal
delta <- truncnorm::rtruncnorm(1, mean = 0, sd = as.numeric(tau[d]), 
    a=0.5 - methPara$mu[changedCpgsIdx] - sqrt(0.25-methPara$var[changedCpgsIdx]), 
    b=0.5 - methPara$mu[changedCpgsIdx] + sqrt(0.25-methPara$var[changedCpgsIdx]))
deltaSim <- c(deltaSim, delta)  # multi core
meaningfulDM <- (abs(delta) >= detectionLimit)
meaningfulDMName <- changedCpgsIdxName[meaningfulDM]

# changing mean methylation for specific CpG's (changedCpgsIdx)
muToBeSimuUNchanged <- methPara$mu[cpgIdx]
muToBeSimuChanged   <- methPara$mu[cpgIdx]
muToBeSimuChanged[match(changedCpgsIdx, cpgIdx)] <- muToBeSimuChanged[match(changedCpgsIdx, cpgIdx)] + delta

## get alpha and beta 
params_unchanged <- getAlphBet(myMean = muToBeSimuUNchanged, myVar = methPara$var[cpgIdx])
alpha_unchanged <- params_unchanged$alpha
beta_unchanged <- params_unchanged$beta
params_changed <- getAlphBet(myMean = muToBeSimuChanged, myVar = methPara$var[cpgIdx])
alpha_changed <- params_changed$alpha
beta_changed <- params_changed$beta

## simulate baseline beta values
g1Beta <- NULL
g2Beta <- NULL
g1Beta <- matrix(stats::rbeta(J*Ncnt, rep(alpha_unchanged, each = Ncnt), rep(beta_unchanged, each = Ncnt)), ncol = Ncnt, byrow = TRUE) 
g2Beta <- matrix(stats::rbeta(J*Ntx, rep(alpha_changed, each = Ntx), rep(beta_changed, each = Ntx)), ncol = Ntx, byrow = TRUE) 
g1Beta[g1Beta == 1] <- max(g1Beta[g1Beta != 1]) # replacing 0/1 by min/max
g2Beta[g2Beta == 1] <- max(g2Beta[g2Beta != 1])
g1Beta[g1Beta == 0] <- min(g1Beta[g1Beta != 0])
g2Beta[g2Beta == 0] <- min(g2Beta[g2Beta != 0])
rownames(g1Beta) <- rownames(g2Beta) <- paste(seq_len(J),"_",names(alpha_unchanged),sep = "")
