#!/usr/bin/env Rscript
# this script does what? MEGENA's example script for plot_module
# except it doesn't work
# well, it was faulty.
library(MEGENA)
library(Cairo)

lv <- ls()

data(Sample_Expression)
# this loads datExpr: 330x844 TGCA subsetted dset
# huge number of samples.

if(length(grep("^ijw$", lv))!=1) {
    ijw <- calculate.correlation(datExpr, doPerm = 2)
}

if(length(grep("^el$", lv))!=1) {
    el <- calculate.PFN(ijw[,1:3])
}
g <- graph.data.frame(el,directed = FALSE)

if(length(grep("^MEGENA.output$", lv))!=1) {
    MEGENA.output <- do.MEGENA(g=g, remove.unsig=F, doPar=F, n.perm = 10)
}


if(length(grep("^output.summary$", lv))!=1) {
output.summary <- MEGENA.ModuleSummary(MEGENA.output,
                                       mod.pvalue = 0.05,hub.pvalue = 0.05,
                                       min.size = 10,max.size = 5000,
                                       annot.table = NULL,id.col = NULL,symbol.col = NULL,
                                       output.sig = TRUE)
}

module.table = output.summary$module.table
colnames(module.table)[1] <- "id"

output.obj <- plot_module_hierarchy(module.table = module.table,
                                    label.scaleFactor=.5,
                                    arrow.size=.025,
                                    edge.color="grey20",
                                    node.label.color="blue")
# print(output.obj[[1]])
CairoPNG("meg3.png", 800, 800)
show(output.obj[[1]])
dev.off()
