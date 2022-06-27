#!/usr/bin/env Rscript
# from 
# https://jdblischak.github.io/dc-bioc-limma/vdx.html
library(Cairo)
library(gridExtra)
library(ggplot2)
library(ggrepel)
library(limma)
library(Biobase)
library(breastCancerVDX)

# Analysis of breast cancer VDX data for videos
# John Blischak
# 2018-08-02
# 
# Study of breast cancer:
# 
#     Bioconductor package breastCancerVDX
#     Published in Wang et al., 2005 and Minn et al., 2007
#     344 patients: 209 ER+, 135 ER-
# 
# Prepare data

data("vdx")
# class(vdx)
# it's an expressionset class.

# dim(vdx)
# Features  Samples 
#    22283      344 

# Phenotypical data:
# pData(vdx)[1:3, ]
#       samplename dataset series id        filename size age er grade pgr
# VDX_3      VDX_3     VDX    VDX  3 GSM36793.CEL.gz   NA  36  0    NA  NA
# VDX_5      VDX_5     VDX    VDX  5 GSM36796.CEL.gz   NA  47  1     3  NA
# VDX_6      VDX_6     VDX    VDX  6 GSM36797.CEL.gz   NA  44  0     3  NA
#       her2 brca.mutation e.dmfs t.dmfs node t.rfs e.rfs treatment tissue
# VDX_3   NA            NA      0   3072    0    NA    NA         0      1
# VDX_5   NA            NA      0   3589    0    NA    NA         0      1
# VDX_6   NA            NA      1    274    0    NA    NA         0      1
#       t.os e.os
# VDX_3   NA   NA
# VDX_5   NA   NA
# VDX_6   NA   NA

#  this is probe data and this is the annotation aspect of it:
# fData(vdx)[1:3, 1:5]
#               probe                                  Gene.title
# 1007_s_at 1007_s_at discoidin domain receptor tyrosine kinase 1
# 1053_at     1053_at replication factor C (activator 1) 2, 40kDa
# 117_at       117_at        heat shock 70kDa protein 6 (HSP70B')
#           Gene.symbol Gene.ID EntrezGene.ID
# 1007_s_at        DDR1     780           780
# 1053_at          RFC2    5982          5982
# 117_at          HSPA6    3310          3310

# short form names of critical variables
x <- exprs(vdx)
f <- fData(vdx)
p <- pData(vdx)

f <- f[, c("Gene.symbol", "EntrezGene.ID", "Chromosome.location")]
colnames(f) <- c("symbol", "entrez", "chrom")

# Recode er as 0 = negative and 1 = positive
p[, "er"] <- ifelse(p[, "er"] == 0, "negative", "positive")
p <- p[, c("id", "age", "er")]

# Explore data
# boxplot(x[1, ] ~ p[, "er"], main = f[1, "symbol"])
eset <- ExpressionSet(assayData = x, phenoData = AnnotatedDataFrame(p), featureData = AnnotatedDataFrame(f))
# dim(eset)
# 
# Features  Samples 
#    22283      344 
# 
# boxplot(exprs(eset)[1, ] ~ pData(eset)[, "er"], main = fData(eset)[1, "symbol"])
# 
# limma pipeline
design <- model.matrix(~er, data = pData(eset))
# head(design, 2)
# 
#       (Intercept) erpositive
# VDX_3           1          0
# VDX_5           1          1
# 
# colSums(design)
# 
# (Intercept)  erpositive 
#         344         209 
# 
# table(pData(eset)[, "er"])
# 
# 
# negative positive 
#      135      209 

fit <- lmFit(eset, design)
head(fit$coefficients, 3)

          (Intercept)  erpositive
1007_s_at   11.725148  0.09878782
1053_at      8.126934 -0.54673000
117_at       7.972049 -0.17342654

fit <- eBayes(fit)
head(fit$t, 3)

          (Intercept) erpositive
1007_s_at    276.8043   1.817824
1053_at      122.5899  -6.428278
117_at       164.0240  -2.781294

results <- decideTests(fit[, "erpositive"])
summary(results)

       erpositive
Down         6276
NotSig      11003
Up           5004

group-means

design <- model.matrix(~0 + er, data = pData(eset))
head(design)

      ernegative erpositive
VDX_3          1          0
VDX_5          0          1
VDX_6          1          0
VDX_7          1          0
VDX_8          1          0
VDX_9          0          1

colSums(design)

ernegative erpositive 
       135        209 

library(limma)
cm <- makeContrasts(status = erpositive - ernegative,
                    levels = design)
cm

            Contrasts
Levels       status
  ernegative     -1
  erpositive      1

fit <- lmFit(eset, design)
head(fit$coefficients)

          ernegative erpositive
1007_s_at  11.725148  11.823936
1053_at     8.126934   7.580204
117_at      7.972049   7.798623
121_at     10.168975  10.086393
1255_g_at   5.903189   5.729195
1294_at     9.166436   9.390949

fit2 <- contrasts.fit(fit, contrasts = cm)
head(fit2$coefficients)

           Contrasts
                 status
  1007_s_at  0.09878782
  1053_at   -0.54673000
  117_at    -0.17342654
  121_at    -0.08258267
  1255_g_at -0.17399402
  1294_at    0.22451339

fit2 <- eBayes(fit2)
results <- decideTests(fit2)
summary(results)

       status
Down     6276
NotSig  11003
Up       5004

topTable(fit2)

              symbol entrez   chrom    logFC   AveExpr        t
205225_at       ESR1   2099  6q25.1 3.762901 11.377735 22.68392
209603_at      GATA3   2625   10p15 3.052348  9.941990 18.98154
209604_s_at    GATA3   2625   10p15 2.431309 13.185334 17.59968
212956_at     TBC1D9  23158 4q31.21 2.157435 11.702942 17.48711
202088_at    SLC39A6  25800 18q12.2 1.719680 13.119496 17.30104
212496_s_at    KDM4B  23030 19p13.3 1.459843 10.703942 16.85070
215867_x_at     CA12    771   15q22 2.246120 11.450485 16.79123
209602_s_at    GATA3   2625   10p15 2.921505  9.547850 16.43202
212195_at      IL6ST   3572    5q11 1.381778 11.737839 16.31864
218195_at   C6orf211  79624  6q25.1 1.738740  9.479901 16.27378
                 P.Value    adj.P.Val         B
205225_at   2.001001e-70 4.458832e-66 149.19866
209603_at   1.486522e-55 1.656209e-51 115.46414
209604_s_at 5.839050e-50 4.337052e-46 102.75707
212956_at   1.665700e-49 9.279201e-46 101.72268
202088_at   9.412084e-49 4.194589e-45 100.01376
212496_s_at 6.188671e-47 2.298369e-43  95.88265
215867_x_at 1.074845e-46 3.421537e-43  95.33780
209602_s_at 3.004184e-45 8.367780e-42  92.05058
212195_at   8.581176e-45 2.124604e-41  91.01458
218195_at   1.299472e-44 2.895613e-41  90.60496

Visualize results

For Ch3 L2 plotMDS/removeBatchEffect

stats <- topTable(fit2, number = nrow(fit2), sort.by = "none")
dim(stats)

[1] 22283     9

hist(runif(10000))

hist(stats[, "P.Value"])

volcanoplot(fit2, highlight = 5, names = fit2$genes[, "symbol"])

Enrichment

topTable(fit2, number = 3)

            symbol entrez  chrom    logFC  AveExpr        t      P.Value
205225_at     ESR1   2099 6q25.1 3.762901 11.37774 22.68392 2.001001e-70
209603_at    GATA3   2625  10p15 3.052348  9.94199 18.98154 1.486522e-55
209604_s_at  GATA3   2625  10p15 2.431309 13.18533 17.59968 5.839050e-50
               adj.P.Val        B
205225_at   4.458832e-66 149.1987
209603_at   1.656209e-51 115.4641
209604_s_at 4.337052e-46 102.7571

# 1000 genes (10% in gene set), 100 are DE (10% in gene set)
fisher.test(matrix(c(10, 100, 90, 900), nrow = 2))


    Fisher's Exact Test for Count Data

data:  matrix(c(10, 100, 90, 900), nrow = 2)
p-value = 1
alternative hypothesis: true odds ratio is not equal to 1
95 percent confidence interval:
 0.4490765 2.0076377
sample estimates:
odds ratio 
         1 

# 1000 genes (10% in gene set), 100 are DE (30% in gene set)
fisher.test(matrix(c(30, 100, 70, 900), nrow = 2))


    Fisher's Exact Test for Count Data

data:  matrix(c(30, 100, 70, 900), nrow = 2)
p-value = 1.88e-07
alternative hypothesis: true odds ratio is not equal to 1
95 percent confidence interval:
 2.306911 6.320992
sample estimates:
odds ratio 
  3.850476 

head(fit2$genes, 3)

          symbol entrez   chrom
1007_s_at   DDR1    780  6p21.3
1053_at     RFC2   5982 7q11.23
117_at     HSPA6   3310    1q23

entrez <- fit2$genes[, "entrez"]


enrich_kegg <- kegga(fit2, geneid = entrez, species = "Hs")
topKEGG(enrich_kegg, number = 4)

                                              Pathway   N Up Down
path:hsa04110                              Cell cycle 115 30   82
path:hsa05166 Human T-cell leukemia virus 1 infection 202 49  124
path:hsa05169            Epstein-Barr virus infection 194 37  114
path:hsa04218                     Cellular senescence 145 32   88
                   P.Up       P.Down
path:hsa04110 0.7955186 8.108269e-16
path:hsa05166 0.9496029 4.456075e-15
path:hsa05169 0.9995466 3.656613e-12
path:hsa04218 0.9784768 9.991751e-11

enrich_go <- goana(fit2, geneid = entrez, species = "Hs")
topGO(enrich_go, ontology = "BP", number = 3)

                            Term Ont    N  Up Down P.Up       P.Down
GO:0002376 immune system process  BP 2337 528 1096    1 1.041713e-41
GO:0006955       immune response  BP 1624 322  797    1 1.403459e-37
GO:0045321  leukocyte activation  BP 1015 216  507    1 1.434672e-25























CairoPNG("pcaplots2.png", 1600, 800)
g0 <- ggplot(xdf1, aes(x=PC1, y=PC2, label=rownames(xdf1))) +geom_point() +geom_text_repel() +
    ggtitle("Batchfx just applied") +
    theme(plot.title = element_text(color="purple", size=14, face="bold.italic"), axis.text=element_text(size=14), axis.title=element_text(size=14,face="bold"))
g1 <- ggplot(xdf2, aes(x=PC1, y=PC2, label=rownames(xdf2))) +geom_point() +geom_text_repel() +
    ggtitle("Batchfx applied and removed") +
    theme(plot.title = element_text(color="purple", size=14, face="bold.italic"), axis.text=element_text(size=14), axis.title=element_text(size=14,face="bold"))
grid.arrange(g0, g1, nrow=1)
dev.off()
