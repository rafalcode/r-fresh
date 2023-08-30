library(SummarizedExperiment)
library(TCGAbiolinks)
library(dplyr)
library(DT) # for the datable function.

# Samples with DNA methylation and gene expression data
# In this example we will access the harmonized database and search for all patients with DNA methylation (platform HumanMethylation450k) and gene expression data for Colon Adenocarcinoma tumor (TCGA-COAD).


call_query <- F
if(call_query) {
    query_meth <- GDCquery(
        project = "TCGA-COAD",
        data.category = "DNA Methylation",
        platform = c("Illumina Human Methylation 450"))
    query_exp <- GDCquery(
        project = "TCGA-COAD",
        data.category = "Transcriptome Profiling",
        data.type = "Gene Expression Quantification", 
        workflow.type = "STAR - Counts")
} else {
    query_exp <- readRDS("coad_query_exp.rds")
    query_meth <- readRDS("coad_query_meth.rds")
}

getres <-  getResults(query_exp, cols = c("data_type","cases"))
expdt <- datatable(getres,
    filter = 'top',
    options = list(scrollX = TRUE, keys = TRUE, pageLength = 5), 
    rownames = FALSE)

brca_met <- GDCquery(
    project = "TCGA-BRCA",
    data.category = "DNA Methylation",
    platform = c("Illumina Human Methylation 450"))

brca_exp <- GDCquery(
    project = "TCGA-BRCA",
    data.category = "Transcriptome Profiling",
    data.type = "Gene Expression Quantification", 
    workflow.type = "STAR - Counts")

# Get all patients that have DNA methylation and gene expression.
common.patients <- intersect(
    substr(getResults(brca_met, cols = "cases"), 1, 12),
    substr(getResults(brca_exp, cols = "cases"), 1, 12))

# Only select the first 5 patients
query_met <- GDCquery(
    project = "TCGA-COAD",
    data.category = "DNA Methylation",
    platform = c("Illumina Human Methylation 450"),
    barcode = common.patients[1:5])

