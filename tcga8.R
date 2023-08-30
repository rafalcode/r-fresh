# Case study n. 1: Pan Cancer downstream analysis BRCA
# frm TCGA biolinks tute numb 8
library(SummarizedExperiment)
library(TCGAbiolinks)

query.exp <- GDCquery(
    project = "TCGA-BRCA", 
    data.category = "Transcriptome Profiling",
    data.type = "Gene Expression Quantification", 
    workflow.type = "STAR - Counts",
    sample.type = c("Primary Tumor","Solid Tissue Normal")
)

GDCdownload(query = query.exp, files.per.chunk = 100)

brca.exp <- GDCprepare(
    query = query.exp, 
    save = TRUE, 
    save.filename = "brcaExp.rda"
)

# get subtype information
infomation.subtype <- TCGAquery_subtype(tumor = "BRCA")

# get clinical data
information.clinical <- GDCquery_clinic(project = "TCGA-BRCA",type = "clinical") 

# Which samples are Primary Tumor
samples.primary.tumour <- brca.exp$barcode[brca.exp$shortLetterCode == "TP"]

# which samples are solid tissue normal
samples.solid.tissue.normal <- brca.exp$barcode[brca.exp$shortLetterCode == "NT"]
Using TCGAnalyze_DEA, we identified 4,815 differentially expression genes (DEG) (log fold change >=1 and FDR < 1%) between 113 normal and 1106 BRCA samples. In order to understand the underlying biological process from DEGs we performed an enrichment analysis using TCGAnalyze_EA_complete function.

dataPrep <- TCGAanalyze_Preprocessing(
    object = brca.exp, 
    cor.cut = 0.6
)                      

dataNorm <- TCGAanalyze_Normalization(
    tabDF = dataPrep,
    geneInfo = geneInfoHT,
    method = "gcContent"
)                

dataFilt <- TCGAanalyze_Filtering(
    tabDF = dataNorm,
    method = "quantile", 
    qnt.cut =  0.25
)   

dataDEGs <- TCGAanalyze_DEA(
    mat1 = dataFilt[,samples.solid.tissue.normal],
    mat2 = dataFilt[,samples.primary.tumour],
    Cond1type = "Normal",
    Cond2type = "Tumor",
    fdr.cut = 0.01 ,
    logFC.cut = 2,
    method = "glmLRT",
    pipeline = "edgeR"
)  
TCGAbiolinks outputs bar chart with the number of genes for the main categories of three ontologies (GO:biological process, GO:cellular component, and GO:molecular function, respectively).

ansEA <- TCGAanalyze_EAcomplete(
    TFname = "DEA genes Normal Vs Tumor",
    RegulonList = dataDEGs$gene_name
)  

TCGAvisualize_EAbarplot(
    tf = rownames(ansEA$ResBP),
    GOBPTab = ansEA$ResBP,
    GOCCTab = ansEA$ResCC,
    GOMFTab = ansEA$ResMF,
    PathTab = ansEA$ResPat,
    nRGTab = dataDEGs$gene_name,
    nBar = 10
)
