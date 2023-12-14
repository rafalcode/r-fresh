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

if(length(grep("^MEGENA.output$", lv))!=1) {
    MEGENA.output <- do.MEGENA(g=g, remove.unsig=F, doPar=F, n.perm = 10)
}

output.summary <- MEGENA.ModuleSummary(MEGENA.output,
                                       mod.pvalue = 0.05,hub.pvalue = 0.05,
                                       min.size = 10,max.size = 5000,
                                       annot.table = NULL,id.col = NULL,symbol.col = NULL,
                                       output.sig = TRUE)

colnames.lst <- c("Saddle Brown","Dark Green","Dark Slate Gray")
# subset module said "comp1_2" I get c1_2 though. ...
pnet.obj <- plot_module(output=output.summary, PFN=g, subset.module="c1_2",
                        # layout="kamada.kawai", label.hubs.only=F,
                        layout="fruchterman.reingold", label.hubs.only=F,
                        gene.set = list("hub.set"=c("CD3E","CD2")), color.code = c("red"),
                        output.plot=F, out.dir="modulePlot",col.names=colnames.lst,
                        hubLabel.col="black",hubLabel.sizeProp=1, show.topn.hubs=Inf,
                        show.legend=T, label.scaleFactor=30, node.sizeProp=23,
                        label.alpha=.8, label.sizeProp=50)
# plot_module(output.summary,PFN,subset.module = NULL,col.names,
# gene.set = NULL,color.code = "logFC",show.legend = TRUE,
# label.hubs.only = TRUE,hubLabel.col = "red",hubLabel.sizeProp = 0.5,show.topn.hubs = 10,
# node.sizeProp = 13,label.sizeProp = 13,label.scaleFactor = 10,label.alpha = 0.5,
# layout = "kamada.kawai",output.plot = TRUE,out.dir = "modulePlot")

# output directory
CairoPNG("pnet.png", 800, 800)
show(pnet.obj)
dev.off()
