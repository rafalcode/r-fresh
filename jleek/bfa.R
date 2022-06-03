# batch effects A chapt from Jeff's Coursera
library(Biobase)
library(sva)
library(bladderbatch)
library(snpStats)

tropical=  c('darkorange', 'dodgerblue', 'hotpink', 'limegreen', 'yellow')
palette(tropical)

data(bladderdata)
pheno = pData(bladderEset)
# "cancer" and "outcome" columns already factors.
# but batch could do with conversion
pheno$batch <- factor(pheno$batch)
edata = exprs(bladderEset)
# edata: num [1:22283, 1:57] 10.12 5.35 6.35 8.9 3.97 ...

## ------------------------------------------------------------------------
modcb = model.matrix(~cancer+batch, data=pheno) # with cancer and batch.
# cancer has three levels, and batch has 5, so giventhere'a base line represented
# by the intercept, we'll have 3-1 + 5-1 +1 for intercept=7 rows in the lm.fit() coefficients matrix.
fit = lm.fit(modcb,t(edata))

# fit is now a list of 8 elements
# "coefficients" "residuals" "effects" "rank" "fitted.values" "assign" "qr" "df.residual"
# view distribution of coefficent for one of the "cancer" levels:
# hist(fit$coefficients[2,],col=2,breaks=100)
# so why does he use row 2 here?
# because it's the "Cancer" level in the cancer factor (yes, I know, confusing)
# note how there will be an intercept in that model, so that's =1 and represents
# the baseline set by the Biopsy level which is fact is the first level for this factor.
# in any case, this is just a graph ... not sure how useful this is, maybe just a sanity check

table(pheno$cancer,pheno$batch)

batch <- pheno$batch
mod0 <- model.matrix(~1, data=pheno) # null model, where you get the sample names coming out as rownames.
# null model, intercept=1 only, "don't touch anything else, don't account for anything else."
modco <- model.matrix(~cancer, data=pheno) # model with cancer var only this time.

# So ComBat() will return the modified edata
combat_edata = ComBat(dat=edata, batch=batch, mod=mod0, par.prior=TRUE, prior.plots=FALSE)
# so we try the lm.fit() again with the batch corrected edata
combat_fit = lm.fit(modcancer,t(combat_edata))
# hist(combat_fit$coefficients[2,],col=2,breaks=100)

# so hs visual inspect of one set of coeeficients, is that they are "smaller"
# he means that a shrinkage has occurred.
# he plots the pre and post combat coefficients against each other:
plot(fit$coefficients[2,],combat_fit$coefficients[2,],col=2,
      xlab="Pre-ComBat LMcoeffs",ylab="Post-ComBat LMcoeffs",xlim=c(-5,5),ylim=c(-5,5))
abline(c(0,1),col=1,lwd=3) # as visual ref .. if there were the same.
# what we get here is a 30degree line more or less showing the new coefficients are smaller
# the intuition for this is that batch effects have been taken out of the matrix
# (all those promoting false positives, to put it in a way) so the effects
# (coefficients) we're seeing are necessarily smaller.

## ------------------------------------------------------------------------
mod = model.matrix(~cancer,data=pheno)
mod0 = model.matrix(~1, data=pheno)
sva1 = sva(edata,mod,mod0,n.sv=2)

## ------------------------------------------------------------------------
summary(lm(sva1$sv ~ pheno$batch))
boxplot(sva1$sv[,2] ~ pheno$batch)
points(sva1$sv[,2] ~ jitter(as.numeric(pheno$batch)),col=as.numeric(pheno$batch))

## ------------------------------------------------------------------------
modsv = cbind(mod,sva1$sv)
fitsv = lm.fit(modsv,t(edata))

## ------------------------------------------------------------------------
par(mfrow=c(1,2))
plot(fitsv$coefficients[2,],combat_fit$coefficients[2,],col=2,
      xlab="SVA",ylab="Combat",xlim=c(-5,5),ylim=c(-5,5))
abline(c(0,1),col=1,lwd=3)
plot(fitsv$coefficients[2,], fit$coefficients[2,],col=2,
      xlab="SVA",ylab="linear model",xlim=c(-5,5),ylim=c(-5,5))
abline(c(0,1),col=1,lwd=3)
