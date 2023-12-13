#!/usr/bin/env Rscript
# this script does what? MEGENA's example script for plot_module
# except it doesn't work
# well, it was faulty.
library(MEGENA)
library(Cairo)

data(Sample_Expression)
# this loads datExpr: 330x844 TGCA subsetted dset
# huge number of samples.

ijw <- calculate.correlation(datExpr[1:100,],doPerm = 2)
el <- calculate.PFN(ijw[,1:3])
g <- graph.data.frame(el,directed = FALSE)

MEGENA.output <- do.MEGENA(g=g, remove.unsig=F, doPar=F, n.perm = 10)

output.summary <- MEGENA.ModuleSummary(MEGENA.output,
                                       mod.pvalue = 0.05,hub.pvalue = 0.05,
                                       min.size = 10,max.size = 5000,
                                       annot.table = NULL,id.col = NULL,symbol.col = NULL,
                                       output.sig = TRUE)

# subset module said "comp1_2" I get c1_2 though. ...
pnet.obj <- plot_module(output=output.summary, PFN=g, subset.module="c1_2",
                        # layout="kamada.kawai", label.hubs.only=F,
                        layout="fruchterman.reingold", label.hubs.only=F,
                        gene.set = list("hub.set"=c("CD3E","CD2")), color.code = c("red"),
                        output.plot=F, out.dir="modulePlot",col.names=c("grey","grey","grey"),
                        hubLabel.col="black",hubLabel.sizeProp=1, show.topn.hubs=Inf,
                        show.legend=T)
# output directory
CairoPNG("pnet.png", 800, 800)
show(pnet.obj)
dev.off()
