#!/usr/bin/env Rscript
# read in that core table gene file which has xls but appears to be a tab delimited thang,
args <- commandArgs(trailingOnly = TRUE)
numargs <- length(args)
enumargs <- 0 # expected numargs
if(numargs != enumargs) {
    print("no args right now")
    stop("Bailing out ..")
}
# ta <- read.table("core_table_gene.xls", sep='\t')
# ta <- read.table("core102.tsv", header=T)
# about 17676A row/genes there.
gq <- grep("pvalue", colnames(ta))
gp <- grep("qvalue", colnames(ta))

# how many genes (rows)
# the compairsons
# "diffexp_deseq2_qvalue_LY2_siCtrl.vs.LY2_siFTO"
# that's col77>
# "diffexp_deseq2_qvalue_LY2_siCtrl.vs.MCF7_siCtrl"
# "diffexp_deseq2_qvalue_LY2_siFTO.vs.MCF7_siFTO"
# "diffexp_deseq2_qvalue_MCF7_siCtrl.vs.MCF7_siFTO"
# "diffexp_deseq2_qvalue_T347_siCtrl.vs.LY2_siCtrl"
# "diffexp_deseq2_qvalue_T347_siCtrl.vs.MCF7_siCtrl"
# "diffexp_deseq2_qvalue_T347_siCtrl.vs.T347_siFTO"
# "diffexp_deseq2_qvalue_T347_siFTO.vs.LY2_siFTO"
# "diffexp_deseq2_qvalue_T347_siFTO.vs.MCF7_siFTO"

# pvalue see how many of each are under the threshold?
lsigp  <- c()
llp <- list() # a list of vector to hold the indices of the sigp genes
cou <- 1
for(i in seq(78,102,3)) {
    llp[[cou]] <- which(ta[i] <0.05)
    lsigp  <- c(lsigp, length(llp[cou]))
    cou <- cou +1
}
names(lsigp) <-c("LY2Ctrl_LY2Fto",
"LY2ctrl_MCF7cTrl",
"LY2fto_MCF7fto",
"MCF7ctrl_MCF7fto",
"T347ctrl_LY2ctrl", 
"T347ctrl_MCF7ctrl", 
"T347ctrl_T347fto",
"T347fto_LY2fto",
"T347fto_MCF7fto")

# qvalue see how many of each are under the threshold?
lsigq  <- c()
llq <- list() # a list of vector to hold the indices of the sigq genes
cou <- 1
for(i in seq(77,101,3)) {
    llq[[cou]] <- which(ta[i] <0.05)
    lsigq  <- c(lsigq, length(llp[cou]))
    cou <- cou +1
}
names(lsigq) <-c("LY2Ctrl_LY2Fto",
"LY2ctrl_MCF7cTrl",
"LY2fto_MCF7fto",
"MCF7ctrl_MCF7fto",
"T347ctrl_LY2ctrl", 
"T347ctrl_MCF7ctrl", 
"T347ctrl_T347fto",
"T347fto_LY2fto",
"T347fto_MCF7fto")

# comments: comparing the ctrls of the three celline 
# gives us large quantities of sig DEGs 10k, and 11k for example
# which is over 50% of all the tested genes.
# the ftos between the 3 is actually similar. Massive DEGs between them.
# Generally when going from p to q there is only a small decline.
#
# what about with-cell-line ctrl vs. fto. Well in MCF7 this is
# very small especially in qvalue, only 88 DEGs here! The p to q effect
# is like losing a factor of 10.
# LY2 aand MCF7 are a bit better, from p to q they lose about 2 thirds of their

