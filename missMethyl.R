# this started actually as a straight wget of the rscript of missMethyl's bioconductor page
## ----load-libs, message=FALSE-------------------------------------------------
library(missMethyl)
library(limma)
library(minfi)
library(Cairo)
library(minfiData)

## ----reading-data, message=FALSE----------------------------------------------
baseDir <- system.file("extdata", package = "minfiData")
targets <- read.metharray.sheet(baseDir)
targets[,1:9]
targets[,10:12]
rgSet <- read.metharray.exp(targets = targets)

## ----ppraw--------------------------------------------------------------------
mSet <- preprocessRaw(rgSet)

## ----swan---------------------------------------------------------------------
mSetSw <- SWAN(mSet,verbose=TRUE)

## ----betasByType, fig.cap = "Beta value dustributions. Density distributions of beta values before and after using SWAN.", echo = TRUE, fig.width=10, fig.height=5----
CairoPNG("betacurves_raw_swan.png", 800, 800)
par(mfrow=c(1,2), cex=1.25)
densityByProbeType(mSet[,1], main = "Raw")
densityByProbeType(mSetSw[,1], main = "SWAN")
dev.off()

## ----filtering----------------------------------------------------------------
detP <- detectionP(rgSet)
keep <- rowSums(detP < 0.01) == ncol(rgSet)
mSetSw <- mSetSw[keep,]

## ----extraction---------------------------------------------------------------
set.seed(10)
mset_reduced <- mSetSw[sample(1:nrow(mSetSw), 20000),]
meth <- getMeth(mset_reduced)
unmeth <- getUnmeth(mset_reduced)
Mval <- log2((meth + 100)/(unmeth + 100))
beta <- getBeta(mset_reduced)
#check: dim(Mval)

## ----mdsplot, fig.cap = "MDS plot. A multi-dimensional scaling (MDS) plot of cancer and normal samples.", echo = TRUE, fig.small=TRUE----
CairoPNG("missMethmds.png", 800, 800)
par(mfrow=c(1,1))
plotMDS(Mval, labels=targets$Sample_Name, col=as.integer(factor(targets$status)))
legend("topleft",legend=c("Cancer","Normal"),pch=16,cex=1.2,col=1:2)
dev.off()

## ----design-------------------------------------------------------------------
group <- factor(targets$status,levels=c("normal","cancer"))
id <- factor(targets$person)
design <- model.matrix(~id + group)
#check: design

## ----diffmeth-----------------------------------------------------------------
fit.reduced <- lmFit(Mval,design)
fit.reduced <- eBayes(fit.reduced, robust=TRUE)

## ----diffmeth-results---------------------------------------------------------
#check: summary(decideTests(fit.reduced))
top<-topTable(fit.reduced,coef=4)
#check: top

## ----top4, fig.cap = "Top DM CpGs. The beta values for the top 4 differentially methylated CpGs.", echo = TRUE, fig.width=10,fig.height=9----
cpgs <- rownames(top)
CairoPNG("missMethtop.png", 800, 800)
par(mfrow=c(2,2))
for(i in 1:4){
    stripchart(beta[rownames(beta)==cpgs[i],]~design[,4],method="jitter",
    group.names=c("Normal","Cancer"),pch=16,cex=1.5,col=c(4,2),ylab="Beta values",
    vertical=TRUE,cex.axis=1.5,cex.lab=1.5)
    title(cpgs[i],cex.main=1.5)
}
dev.off()

## ----diffmeth2----------------------------------------------------------------
# get M-values for ALL probes
meth <- getMeth(mSet)
unmeth <- getUnmeth(mSet)
M <- log2((meth + 100)/(unmeth + 100))
# setup the factor of interest
grp <- factor(targets$status, labels=c(0,1))
# extract Illumina negative control data
INCs <- getINCs(rgSet)
#check: head(INCs)
# add negative control data to M-values
Mc <- rbind(M,INCs)
# create vector marking negative controls in data matrix
ctl1 <- rownames(Mc) %in% rownames(INCs)
#check: table(ctl1)
rfit1 <- RUVfit(Y = Mc, X = grp, ctl = ctl1) # Stage 1 analysis
rfit2 <- RUVadj(Y = Mc, fit = rfit1)

## ----ruv1---------------------------------------------------------------------
top1 <- topRUV(rfit2, num=Inf, p.BH = 1)
#check: head(top1)
ctl2 <- rownames(M) %in% rownames(top1[top1$p.BH_X1.1 > 0.5,])
#check: table(ctl2)

## ----ruv2---------------------------------------------------------------------
# Perform RUV adjustment and fit
rfit3 <- RUVfit(Y = M, X = grp, ctl = ctl2) # Stage 2 analysis
rfit4 <- RUVadj(Y = M, fit = rfit3)
# Look at table of top results
#check: topRUV(rfit4)

## ----limmaruv-----------------------------------------------------------------
# setup design matrix
des <- model.matrix(~grp)
#check: des
# limma differential methylation analysis
lfit1 <- lmFit(M, design=des)
lfit2 <- eBayes(lfit1) # Stage 1 analysis
# Look at table of top results
#check: topTable(lfit2)

## ----limmaruv1----------------------------------------------------------------
topl1 <- topTable(lfit2, num=Inf)
#check: head(topl1)
ctl3 <- rownames(M) %in% rownames(topl1[topl1$adj.P.Val > 0.5,])
#check: table(ctl3)

## ----limmaruv2----------------------------------------------------------------
# Perform RUV adjustment and fit
rfit5 <- RUVfit(Y = M, X = grp, ctl = ctl3) # Stage 2 analysis
rfit6 <- RUVadj(Y = M, fit = rfit5)
# Look at table of top results
#check: topRUV(rfit6)

## ----ruvadj-------------------------------------------------------------------
Madj <- getAdj(M, rfit5) # get adjusted values

## ----mdsplotadj, fig.cap = "RUVm adjusted data. An MDS plot of cancer and normal data, before and after RUVm adjustment.", echo = TRUE, fig.width=10, fig.height=5----
CairoPNG("missMethmdsplotadj.png", 800, 800)
par(mfrow=c(1,2))
plotMDS(M, labels=targets$Sample_Name, col=as.integer(factor(targets$status)),
        main="Unadjusted", gene.selection = "common")
legend("right",legend=c("Cancer","Normal"),pch=16,cex=1,col=1:2)
plotMDS(Madj, labels=targets$Sample_Name, col=as.integer(factor(targets$status)),
        main="Adjusted: RUV-inverse", gene.selection = "common")
legend("topright",legend=c("Cancer","Normal"),pch=16,cex=1,col=1:2)
dev.off()

## ----ruvadj1------------------------------------------------------------------
# Use RUV-4 in stage 2 of RUVm with k=1 and k=2
rfit7 <- RUVfit(Y = M, X = grp, ctl = ctl3,
                method = "ruv4", k=1) # Stage 2 with RUV-4, k=1
rfit9 <- RUVfit(Y = M, X = grp, ctl = ctl3,
                method = "ruv4", k=2) # Stage 2 with RUV-4, k=2
# get adjusted values
Madj1 <- getAdj(M, rfit7)
Madj2 <- getAdj(M, rfit9)

## ----mdsplotadj1, fig.cap = "Effect of different adjustment methods and parameters. MDS plots of cancer and normal data before an after adjustment with RUV-inverse and RUV-4 with different k values.", echo = TRUE, fig.width=10, fig.height=9----
CairoPNG("missMethmdsplotadj1.png", 800, 800)
par(mfrow=c(2,2))
plotMDS(M, labels=targets$Sample_Name, col=as.integer(factor(targets$status)),
        main="Unadjusted", gene.selection = "common")
legend("top",legend=c("Cancer","Normal"),pch=16,cex=1,col=1:2)
plotMDS(Madj, labels=targets$Sample_Name, col=as.integer(factor(targets$status)),
        main="Adjusted: RUV-inverse", gene.selection = "common")
legend("topright",legend=c("Cancer","Normal"),pch=16,cex=1,col=1:2)
plotMDS(Madj1, labels=targets$Sample_Name, col=as.integer(factor(targets$status)),
        main="Adjusted: RUV-4, k=1", gene.selection = "common")
legend("bottom",legend=c("Cancer","Normal"),pch=16,cex=1,col=1:2)
plotMDS(Madj2, labels=targets$Sample_Name, col=as.integer(factor(targets$status)),
        main="Adjusted: RUV-4, k=2", gene.selection = "common")
legend("bottomright",legend=c("Cancer","Normal"),pch=16,cex=1,col=1:2)
dev.off()

## ----checkdesign: design

## ----diffvar------------------------------------------------------------------
fitvar <- varFit(Mval, design = design, coef = c(1,4))

## ----diffvar-results----------------------------------------------------------
#check: summary(decideTests(fitvar))
topDV <- topVar(fitvar, coef=4)
#check:topDV

## ----alternative--------------------------------------------------------------
design2 <- model.matrix(~0+group+id)
fitvar.contr <- varFit(Mval, design=design2, coef=c(1,2))
contr <- makeContrasts(groupcancer-groupnormal,levels=colnames(design2))
fitvar.contr <- contrasts.varFit(fitvar.contr,contrasts=contr)

## ----altresults---------------------------------------------------------------
#check: summary(decideTests(fitvar.contr))
#check: topVar(fitvar.contr,coef=1)

## ----top4DV,fig.cap="Top DV CpGs. The beta values for the top 4 differentially variable CpGs.", fig.width=10, fig.height=9----
cpgsDV <- rownames(topDV)
CairoPNG("missMethtopdv.png", 800, 800)
par(mfrow=c(2,2))
for(i in 1:4){
stripchart(beta[rownames(beta)==cpgsDV[i],]~design[,4],method="jitter",
group.names=c("Normal","Cancer"),pch=16,cex=1.5,col=c(4,2),ylab="Beta values",
vertical=TRUE,cex.axis=1.5,cex.lab=1.5)
title(cpgsDV[i],cex.main=1.5)
}
dev.off()

## ----loadingdata--------------------------------------------------------------
library(tweeDEseqCountData)
data(pickrell1)
counts<-exprs(pickrell1.eset)
#check: dim(counts)
gender <- pickrell1.eset$gender
#check: table(gender)
rm(pickrell1.eset)
data(genderGenes)
data(annotEnsembl63)
annot <- annotEnsembl63[,c("Symbol","Chr")]
rm(annotEnsembl63)

## ----dgelist------------------------------------------------------------------
library(edgeR)
y <- DGEList(counts=counts, genes=annot[rownames(counts),])

## ----dgelist-filtering--------------------------------------------------------
isexpr <- rowSums(cpm(y)>1) >= 20
hasannot <- rowSums(is.na(y$genes))==0
y <- y[isexpr & hasannot,,keep.lib.sizes=FALSE]
#check: dim(y)
y <- calcNormFactors(y)

## ----testhapmap---------------------------------------------------------------
design.hapmap <- model.matrix(~gender)
fitvar.hapmap <- varFit(y, design = design.hapmap, coef=c(1,2))
fitvar.hapmap$genes <- y$genes

## ----resultshapmap------------------------------------------------------------
#check: summary(decideTests(fitvar.hapmap))
topDV.hapmap <- topVar(fitvar.hapmap,coef=ncol(design.hapmap))
#check: topDV.hapmap

## ----top4DVhapmap,fig.cap="Top DV CpGs. The log counts per million for the top 4 differentially variably expressed genes.", fig.width=10, fig.height=9----
genesDV <- rownames(topDV.hapmap)
CairoPNG("missMethtop4dvhap.png", 800, 800)
par(mfrow=c(2,2))
for(i in 1:4){
    stripchart(cpm(y,log=TRUE)[rownames(y)==genesDV[i],]~design.hapmap[,ncol(design.hapmap)],method="jitter",
    group.names=c("Female","Male"),pch=16,cex=1.5,col=c(4,2),ylab="Log counts per million",
    vertical=TRUE,cex.axis=1.5,cex.lab=1.5)
    title(genesDV[i],cex.main=1.5)
}
dev.off()

## ----gometh1------------------------------------------------------------------
top <- topRUV(rfit4, number = Inf, p.BH = 1)
#check: table(top$p.BH_X1.1 < 0.01)

## ----gometh2------------------------------------------------------------------
beta <- getBeta(mSet)
# make sure that order of beta values matches orer after analysis
beta <- beta[match(rownames(top),rownames(beta)),]
beta_norm <- rowMeans(beta[,grp==0])
beta_can <- rowMeans(beta[,grp==1])
Delta_beta <- beta_can - beta_norm
sigDM <- top$p.BH_X1.1 < 0.01 & abs(Delta_beta) > 0.25
#check: table(sigDM)

## ----gometh3------------------------------------------------------------------
topCpGs<-topRUV(rfit4,number=10000)
sigCpGs <- rownames(topCpGs)
#check: sigCpGs[1:10]

# Check number of genes that significant CpGs are annotated to
check <- getMappedEntrezIDs(sig.cpg = sigCpGs)
#check: length(check$sig.eg)

## ----gometh4, fig.cap="Probe number bias in the cancer dataset.", fig.width=6, fig.height=5----
library(IlluminaHumanMethylation450kanno.ilmn12.hg19)
gst <- gometh(sig.cpg=sigCpGs, all.cpg=rownames(top), collection="GO", plot.bias=F)
#check: topGSA(gst, n=10)

## ----gometh5------------------------------------------------------------------
gst.kegg <- gometh(sig.cpg=sigCpGs, all.cpg=rownames(top), collection="KEGG")
#check: topGSA(gst.kegg, n=10)

## ----gometh6------------------------------------------------------------------
gst.kegg.prom <- gometh(sig.cpg=sigCpGs, all.cpg=rownames(top), 
                        collection="KEGG", genomic.features = c("TSS200",
                                                                "TSS1500",
                                                                "1stExon"))
#check: topGSA(gst.kegg.prom, n=10)

## ----gometh7------------------------------------------------------------------
gst.kegg.body <- gometh(sig.cpg=sigCpGs, all.cpg=rownames(top), 
                        collection="KEGG", genomic.features = c("Body"))
#check: topGSA(gst.kegg.body, n=10)

## ----gometh8------------------------------------------------------------------
gst.kegg.body <- gometh(sig.cpg=sigCpGs, all.cpg=rownames(top), 
                        collection="KEGG", genomic.features = c("Body"), 
                        sig.genes = TRUE)
#check: topGSA(gst.kegg.body, n=5)

## ----gsameth------------------------------------------------------------------
hallmark <- readRDS(url("http://bioinf.wehi.edu.au/MSigDB/v7.1/Hs.h.all.v7.1.entrez.rds"))
gsa <- gsameth(sig.cpg=sigCpGs, all.cpg=rownames(top), collection=hallmark)
#check: topGSA(gsa, n=10)

## ----dmrcate1-----------------------------------------------------------------
library(DMRcate)

## ----dmrcate2-----------------------------------------------------------------
myAnnotation <- cpg.annotate(object = M, datatype = "array", what = "M", 
                             arraytype = c("450K"), 
                             analysis.type = "differential", design = design, 
                             coef = 4)

## ----dmrcate3-----------------------------------------------------------------
DMRs <- dmrcate(myAnnotation, lambda=1000, C=2)
results.ranges <- extractRanges(DMRs)
#check: results.ranges

## ----dmrcatetopDMR,fig.cap="Top DMR from DMRcate.", fig.width=10, fig.height=9----
cols <- c(2,4)[group]
names(cols) <-group
beta <- getBeta(mSet)

gettinginternalservererror <- T
if(!gettinginternalservererror) {
    CairoPNG("missMethdmrcatetop.png", 800, 800)
    par(mfrow=c(1,1))
    DMR.plot(ranges=results.ranges, dmr=2, CpGs=beta, phen.col=cols, 
             what="Beta", arraytype="450K", genome="hg19")
    dev.off()
} 

## ----goregion1, fig.cap="Probe number bias for DMRs in the cancer dataset.", fig.width=6, fig.height=5----
gst.region <- goregion(results.ranges, all.cpg=rownames(M), 
                       collection="GO", array.type="450K", plot.bias=F)

## ----goregion2----------------------------------------------------------------
#check: topGSA(gst.region, n=10)

## ----goregion3----------------------------------------------------------------
gst.region.kegg <- goregion(results.ranges, all.cpg=rownames(M), 
                       collection="KEGG", array.type="450K")
#check: topGSA(gst.region.kegg, n=10)

## ----gsaregion----------------------------------------------------------------
gsa.region <- gsaregion(results.ranges, all.cpg=rownames(M), 
                        collection=hallmark)
#check: topGSA(gsa.region, n=10)
#FINIS
