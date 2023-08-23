#!/usr/bin/env Rscript
# this script does what? resex3.R, the example 3 from the helpfile of DESeq2's results() function.
# I've always found DESeq2's vignette a disappointing read.
# Don't get me wrong, it goes to great effort,
# but it semes to ger cuaght up in itself.
# the section onf interactions is fluffy .. is not incisive.
# in any case it actually refers you to the manpage for the hideously named results() function.
# that is replicated here.
library(DESeq2)
library(Cairo)

## Example 3: two conditions, three genotypes
# subtitle (I guess)
# ~~~ Using interaction terms ~~~

dds <- makeExampleDESeqDataSet(n=100,m=18) # this already has two conditions built in , first 9 are A, second are B.

dds$genotype <- factor(rep(rep(c("I","II","III"),each=3),2))

# design(dds) <- ~ genotype + condition + genotype:condition
design(dds) <- ~ condition*genotype

dds <- DESeq(dds)
# check:
resna <- resultsNames(dds)
# you get
# [1] "Intercept"              "genotype_II_vs_I"       "genotype_III_vs_I"      "condition_B_vs_A"       "genotypeII.conditionB"
# [6] "genotypeIII.conditionB" 

# the condition effect for genotype I (the main effect)
# thi sis the bog standard, one first factor (main effect is considered)
# it doesn't really pertain to the specific problem we have here.
# so it's a bit of a throwaway line.
res0 <- results(dds, contrast=c("condition","B","A"))

# the condition effect for genotype III.
# this is the main effect *plus* the interaction term
# (the extra condition effect in genotype III compared to genotype I).
# res1 <- results(dds, contrast=list( c("condition_B_vs_A","genotypeIII.conditionB") ))
# I wonder if this will work:
res1 <- results(dds, contrast=list(resna[2], resna[6])) # yes, it does.

# the interaction term for condition effect in genotype III vs genotype I.
# this tests if the condition effect is different in III compared to I
res2 <- results(dds, name="genotypeIII.conditionB")

# the interaction term for condition effect in genotype III vs genotype II.
# this tests if the condition effect is different in III compared to II
res4 <- results(dds, contrast=list("genotypeIII.conditionB", "genotypeII.conditionB"))

# Note that a likelihood ratio could be used to test if there are any
# differences in the condition effect between the three genotypes.

# ~~~ Using a grouping variable ~~~

# This is a useful construction when users just want to compare
# specific groups which are combinations of variables.
dds2 <- dds
dds2$group <- factor(paste0(dds$genotype, dds$condition))

design(dds2) <- ~ group
dds2 <- DESeq(dds2)
resultsNames(dds2)

# the condition effect for genotypeIII
res5 <- results(dds2, contrast=c("group", "IIIB", "IIIA"))

# Cairo image template
# CairoPNG("fname.png", 800, 800)
# put plot command here
# dev.off()
