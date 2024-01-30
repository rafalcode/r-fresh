#!/usr/bin/env Rscript
# bmrt0.R .. reflects biomaRt ops
# this script does what? Various biomaRt tests
# especially where versions are concerned
# and the infanmous "no httr slot" error.
# which has to do with old Mart objects not having that slot, and the new biomaRt not liking that.
library(biomaRt)

# ensembl <- useEnsembl(biomart = "genes", dataset = "hsapiens_gene_ensembl")
# this actually gives a Mart object
# with 8 slots, the last being httr_config

# ensembl <- useEnsembl(biomart = "genes", dataset = "hsapiens_gene_ensembl", host="http://dec2017.archive.ensembl.org")
# there appears to be no dec2017 server anymore.
# ensembl <- useEnsembl(biomart = "genes", dataset = "hsapiens_gene_ensembl", host="https://jul2018.archive.ensembl.org")
ensembl <- useEnsembl(biomart = "genes", dataset = "hsapiens_gene_ensembl", host="https://may2015.archive.ensembl.org")

# FUNNILY enough, the may2015 DOES have the httr_config slot.
