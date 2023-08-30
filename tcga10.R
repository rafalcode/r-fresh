# fromt TCGA biolonks tute numb 10
# Classifiers methods
# oh hang on this is for glioma only 
# a whole function devoted to calissifyng them
# fine but opaque
# 
# 
# Classifying gliomas samples with gliomaClassifier
# Classifying glioma samples with DNA methylation array based on:
# 
# Ceccarelli, Michele, et al. “Molecular profiling reveals biologically discrete subsets and pathways of progression in diffuse glioma.” Cell 164.3 (2016): 550-563. (https://doi.org/10.1016/j.cell.2015.12.028)
# 
# Possible classifications are:
# 
# Mesenchymal-like
# Classic-like
# G-CIMP-high
# G-CIMP-low
# LGm6-GBM
# Codel
# Data
# The input data can be either a Summarized Experiment object of a matrix (samples as columns, probes as rows) from the following platforms:
# 
# HM27
# HM450
# EPIC array.
# In this example we will retrieve two samples from TCGA and classify them expecting the same result as the paper.

library(SummarizedExperiment)
library(TCGAbiolinks)

query <- GDCquery(
  project = "TCGA-GBM",
  data.category = "DNA methylation",
  barcode = c("TCGA-06-0122","TCGA-14-1456"),
  platform = "Illumina Human Methylation 27",
  legacy = TRUE
)

GDCdownload(query)
data.hg19 <- GDCprepare(query)
# assay(data.hg19)[1:5,1:2]

# Function
classification <- gliomaClassifier(data.hg19)

# Results
# The classfier will return a list of 3 data frames:

Sample final classification
Each model final classification
Each class probability of classification
names(classification)
classification$final.classification
classification$model.classifications
classification$model.probabilities
Comparing results with paper
TCGAquery_subtype("GBM") %>%
 dplyr::filter(patient %in% c("TCGA-06-0122","TCGA-14-1456")) %>%
 dplyr::select("patient","Supervised.DNA.Methylation.Cluster")
## gbm subtype information from:doi:10.1016/j.cell.2015.12.028
patient
<fct>
Supervised.DNA.Methylation.Cluster
<fct>
TCGA-06-0122	Classic-like
TCGA-14-1456	G-CIMP-high
2 rows
