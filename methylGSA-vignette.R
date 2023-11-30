library(methylGSA)
library(IlluminaHumanMethylation450kanno.ilmn12.hg19)
library(org.Hs.eg.db)
library(BiocParallel)

data(cpgtoy)

## -----------------------------------------------------------------------------
res1 = methylglm(cpg.pval = cpg.pval, minsize = 200, maxsize = 500, GS.type = "KEGG")
#check: head(res1, 15)

glm_res = data.frame(
    "Column" = c("ID", 
                 "Description (N/A for user supplied gene set)",
                 "Size", "pvalue", "padj"),
    "Explanation" = c("Gene set ID", "Gene set description", 
                "Number of genes in gene set", 
                "p-value in logistic regression", "Adjusted p-value")
)

## -----------------------------------------------------------------------------
genes_04080 = select(org.Hs.eg.db, "04080", "SYMBOL", keytype = "PATH")
#check head(genes_04080)

#  # include all the IDs as the 2nd argument in select function
genes_all_pathway = select(org.Hs.eg.db, as.character(res1$ID), "SYMBOL", keytype = "PATH")
#check  head(genes_all_pathway)

res2 <- methylRRA(cpg.pval = cpg.pval, method = "ORA", minsize = 200, maxsize = 210)
#check:  head(res2, 15)

## ----echo=FALSE---------------------------------------------------------------
ora_res <- data.frame(
    "Column" = c("ID", 
                 "Description (N/A for user supplied gene set)", "Count",
                 "overlap", "Size", "pvalue", "padj"),
    "Explanation" = c("Gene set ID", "Gene set description", 
                "Number of significant genes in the gene set",
                "Names of significant genes in the gene set",
                "Number of genes in gene set", 
                "p-value in ORA", "Adjusted p-value")
)

res3 <- methylRRA(cpg.pval = cpg.pval, method = "GSEA", minsize = 200, maxsize = 210)

gsea_res = data.frame(
    "Column" = c("ID", 
                 "Description (N/A for user supplied gene set)", 
                 "Size", "enrichmentScore", "NES",
                 "pvalue", "padj"),
    "Explanation" = c("Gene set ID", "Gene set description", 
                "Number of genes in gene set",
                "Enrichment score (see [3] for details)", 
                "Normalized enrichment score (see [3] for details)",
                "p-value in GSEA", "Adjusted p-value")
)

res4 = methylgometh(cpg.pval = cpg.pval, sig.cut = 0.001, minsize = 200, maxsize = 210)

data(GSlisttoy)
## to make the display compact, only a proportion of each gene set is shown
# head(lapply(GS.list, function(x) x[1:30]), 3)   

res_p = methylglm(cpg.pval = cpg.pval, minsize = 200, maxsize = 500, GS.type = "KEGG", parallel = TRUE)

## -----------------------------------------------------------------------------
data(CpG2Genetoy)
# head(CpG2Gene)   
FullAnnot = prepareAnnot(CpG2Gene) 

## -----------------------------------------------------------------------------
GS.list = GS.list[1:10]
res5 = methylRRA(cpg.pval = cpg.pval, FullAnnot = FullAnnot, method = "ORA", 
                    GS.list = GS.list, GS.idtype = "SYMBOL", 
                    minsize = 100, maxsize = 300)

res6 = methylglm(cpg.pval = cpg.pval, array.type = "450K", GS.type = "Reactome", minsize = 100, maxsize = 110)

## -----------------------------------------------------------------------------
CairoPNG("mgsabar.png", 800, 800)
barplot(res1, num = 8, colorby = "pvalue")
dev.off()
