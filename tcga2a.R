library(SummarizedExperiment)
library(TCGAbiolinks)
library(dplyr)
library(DT) # for the datable function.

# Samples with DNA methylation and gene expression data
# In this example we will access the harmonized database and search for all patients with DNA methylation (platform HumanMethylation450k) and gene expression data for Colon Adenocarcinoma tumor (TCGA-COAD).


brca_exp <- GDCquery(
    project = "TCGA-BRCA",
    data.category = "Transcriptome Profiling",
    data.type = "Gene Expression Quantification", 
    workflow.type = "STAR - Counts")

# Get all patients that have DNA methylation and gene expression.
getcases <- getResults(brca_exp, cols = "cases")

GDCdownload(brca_exp)
data <- GDCprepare(brca_exp)
