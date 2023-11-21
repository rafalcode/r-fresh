# you're right this is about combat
# remember
# only one batch ( Ihtnk this means batch variable, only 1 batch vale doesn't make sense, there's got to be more than one).
# 
# however, despite the name this ended up not being abou tcombat at all.
# comba.R was getting complicated.
library(Cairo)
library(ggrepel) # this brings in ggplot
library(ggforce)
library(sva)
library(bladderbatch)
library(stargazer)
library(jtools)
library(Polychrome)
library(limma)

plotPCArflev3 <- function(ldf, pngname="mypcaplot.png")
{
    # allow not only for extra factor, but also cuttree output
    # note: no pca calc here .. it must be precalculated
    require(Cairo)
    require(ggrepel)
    require(Polychrome)

    # factor fac, must be linked to mat, same number of columns
    # yes, now they are in the df this is guaranteed.

    #colours
    levsfac <- nlevels(ldf[[1]]$fac)
    firstpal <- Polychrome::createPalette(levsfac, c("#010101", "#ff0000"), M=1000)

    CairoPNG(pngname, 1000, 800)
        # ggp <- ggplot(data=df, aes(x=!!sym("PC1"), y=!!sym("PC2"), col=!!sym("fac"), label=!!sym("name"))) +
        ggp <- ggplot(data=ldf[[1]], aes(x=!!sym("PC1"), y=!!sym("PC2"), col=!!sym("fac"))) +
        xlim(-120,100) +
        geom_point(size=3) +
        xlab(ldf[[2]]) + ylab(ldf[[3]]) + 
        geom_text_repel(label=ldf[[1]]$name, size=3, position = position_nudge_repel(x = 0.1, y = 0.1)) +
        coord_fixed() +
        theme_minimal() +
        # ggforce::geom_mark_hull(aes(group = cutt), colour="black")
        # ggforce::geom_mark_hull(aes(group=df$cutt, label=df$cutt), colour="grey", concavity=.25, show.legend=T)
        ggforce::geom_mark_hull(aes(group=ldf[[1]]$cutt, label=ldf[[1]]$cutt), colour="grey", concavity=.25)
    show(ggp)
    dev.off()
} # end of funct

plotPCArflev2 <- function(df, fac, cutt, pngname="mypcaplot.png")
{
    # allow not only for extra factor, but also cuttree output
    # note: no pca calc here .. it must be precalculated
    require(Cairo)
    require(ggrepel)
    require(Polychrome)
    # factor fac, must be linked to mat, same number of columns
    lfac <- length(fac)
    if(lfac != nrow(df)) {
        stop("Don't use this function if length(fac) != ncol(mat)\n")
    }
    lcutt <- length(cutt)
    if(lcutt != nrow(df)) {
        stop("Don't use this function if length(cuttree) != ncol(mat)\n")
    }
    #colours
    levsfac <- nlevels(fac)
    firstpal <- Polychrome::createPalette(levsfac, c("#010101", "#ff0000"), M=1000)

    CairoPNG(pngname, 1000, 800)
    ggp <- ggplot(data=df, aes(x=!!sym("PC1"), y=!!sym("PC2"), col=fac, label=!!sym("name"))) +
        geom_point(size=3) +
        xlab(paste0("PC1: ",round(df$percentVar[1] * 100),"% variance")) +
        ylab(paste0("PC2: ",round(df$percentVar[2] * 100),"% variance")) +
        geom_text_repel(size=4, position = position_nudge_repel(x = 0.1, y = 0.1)) +
        coord_fixed() +
        theme_minimal()
    show(ggp)
    dev.off()
} # end of funct

plotPCArflev <- function(mat, fac, pngname="mypcaplot.png")
{
    # allow for extra factor that matches columns
    require(Cairo)
    require(ggrepel)
    require(Polychrome)
    # factor fac, must be linked to mat, same number of columns
    lfac <- length(fac)
    if(lfac != ncol(mat)) {
        stop("Don;t use this function if length(fac) != ncol(mat)\n")
    }
    #colours
    levsfac <- nlevels(fac)
    firstpal <- Polychrome::createPalette(levsfac, c("#010101", "#ff0000"), M=1000)
    names(firstpal) <- NULL
    pca <- prcomp(t(mat)) # cuz we are expecting a count matrix

    # the contribution to the total variance for each component
    percentVar <- pca$sdev^2 / sum( pca$sdev^2 )
    d <- data.frame(PC1=pca$x[,1], PC2=pca$x[,2], name=colnames(mat))

    CairoPNG(pngname, 800, 800)
    ggp <- ggplot(data=d, aes(x=!!sym("PC1"), y=!!sym("PC2"), col=fac, label=!!sym("name"))) +
        geom_point(size=3) +
        xlab(paste0("PC1: ",round(percentVar[1] * 100),"% variance")) +
        ylab(paste0("PC2: ",round(percentVar[2] * 100),"% variance")) +
        geom_text_repel(size=4, position = position_nudge_repel(x = 0.1, y = 0.1)) +
        # coord_fixed()
        # meant to work but doesn't (try stripping colours from the colourhex vector:
        coord_fixed() +
        scale_color_manual(values=firstpal)
    show(ggp)
    dev.off()
    return(d)
} # end of funct

couma2pca <- function(coumat)
{
    pca <- prcomp(t(coumat)) # cuz we are expecting a count matrix

    # the contribution to the total variance for each component
    percentVar <- pca$sdev^2 / sum( pca$sdev^2 )
    d <- data.frame(PC1=pca$x[,1], PC2=pca$x[,2], name=colnames(coumat))
    pv1lab <- xlab(paste0("PC1:",round(percentVar[1] * 100),"%var"))
    pv2lab <- xlab(paste0("PC2:",round(percentVar[2] * 100),"%var"))
    return(list(d, pv1lab, pv2lab))
} # end of funct

lv <- ls()
data(bladderdata)

pheno <- pData(bladderEset)
pheno$batch <- factor(pheno$batch)

# The expression data can be obtained from the expression slot of the expression set.
edata <- exprs(bladderEset)
colnames(edata) <- gsub("^GSM710([^\\.]+)\\.CEL$", "\\1",colnames(edata))

# d <- plotPCArflev(edata, pheno$cancer, "eda1.png")
ldf <- couma2pca(edata) # we're getting a list here, a df and two labels.

dd <- dist(ldf[[1]])
hc <- hclust(dd)
CairoPNG("bb0hc.png", 800, 800)
plot(hc)
dev.off()

cutt0 <- cutree(hc, 6)
names(cutt0) <- NULL
cutt <- paste0("G", as.character(cutt0))
ldf[[1]]$cutt <- factor(cutt)
ldf[[1]]$fac <- pheno$cancer
plotPCArflev3(ldf, "eda2.png")

# mod <- model.matrix(~as.factor(cancer), data=pheno)
# as.factor() not required, already is one.
# mm <- model.matrix(~ cancer + outcome, data=pheno) # manually you can see outcome closely related to cancer var.
# in any case it prodcues singualrities
mm <- model.matrix(~ cancer,  data=pheno)

# The null model contains only the adjustment variables. Since we are not adjusting for any other variables in this analysis, only an intercept is included in the model.
m0 <- model.matrix(~1,data=pheno)

if(1!=length(grep("^svafit$", lv))) {
    svafit <- sva(edata, mod=mm, mod0=m0)
    numsv <- ncol(svafit$sv)
}

iwantosee <- T
if(iwantosee) {
    for (i in 1:numsv) {
        y  <- svafit$sv[, i]
        svlabel <- paste0("SV", i)

        m1 <- lm(y ~ 0 + batch, data = pheno)
        m2 <- lm(y ~ 0 + cancer, data = pheno)
        m3 <- lm(y ~ 0 + outcome , data = pheno)
        m4 <- lm(y ~ batch + cancer + outcome, data = pheno)
        col_labs <- c("Batch", "Cancer", "Outcome", "All")

        # output the linear model results for surrogate variables and batch
        stargazer::stargazer(m1, m2, m3, m4, dep.var.labels = svlabel,
                             column.labels = col_labs, single.row = TRUE, type = "text", title = "Batch Effect Assessment",
                             out = paste0("bladderbatch-sv-linear-model-results-", svlabel, ".txt"))

        CairoPNG(paste0("bladdb.sv", i, ".psumms.png"), 800, 800)
        pj <- jtools::plot_summs(m1, m2, m3, m4, model.names = col_labs, legend.title = svlabel)
        show(pj)
        dev.off()
    }
}


edat2 <- limma::removeBatchEffect(edata, covariates = svafit$sv[,1])
ld2 <- couma2pca(edat2) # we're getting a list here, a df and two labels.

d2 <- dist(ld2[[1]])
hc2 <- hclust(d2)
CairoPNG("bb2hc.png", 800, 800)
plot(hc2)
dev.off()

cutt0 <- cutree(hc2, 6)
names(cutt0) <- NULL
cutt <- paste0("G", as.character(cutt0))
ld2[[1]]$cutt <- factor(cutt)
ld2[[1]]$fac <- pheno$cancer
plotPCArflev3(ld2, "eda3.png")
