library(Cairo)

CairoPDF("ircplot.pdf", 6, 6, bg="transparent")
data(iris)
attach(iris)
plot(Petal.Length, rep(-0.03,length(Species)), xlim=c(1,7),
     ylim=c(0,1.7), xlab="Petal.Length", ylab="Density",
     pch=21, cex=1.5, col="#00000001", main = "Iris (yet again)",
     bg=c("#ff000020","#00ff0020","#0000ff20")[unclass(Species)])
for (i in 1:3)
  polygon(density(Petal.Length[unclass(Species)==i],bw=0.2),
    col=c("#ff000040","#00ff0040","#0000ff40")[i])
dev.off()

pdf("irpplot.pdf", bg="transparent", paper="A4")
data(iris)
attach(iris)
plot(Petal.Length, rep(-0.03,length(Species)), xlim=c(1,7),
     ylim=c(0,1.7), xlab="Petal.Length", ylab="Density",
     pch=21, cex=1.5, col="#00000001", main = "Iris (yet again)",
     bg=c("#ff000020","#00ff0020","#0000ff20")[unclass(Species)])
for (i in 1:3)
  polygon(density(Petal.Length[unclass(Species)==i],bw=0.2),
    col=c("#ff000040","#00ff0040","#0000ff40")[i])
dev.off()
