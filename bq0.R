#!/usr/bin/env Rscript
# this script does what?
# from https://stackoverflow.com/questions/58048861/legend-using-bquote-together-with-a-vector-variable
library(Cairo)

kappa_var <- c(0.5, 1.0, 1.5, 5.0)
kappa_var <- c("0.5", "1.0", "1.5", "5.0")

CairoPNG("bq0.png", 800, 800)

plot(1000, type="n", xlab="x", ylab= expression(paste("f(x;",kappa,",",sigma,")",sep="")), 
     xlim=c(0, 5), ylim=c(0, 6))

# legend("topleft",legend=do.call( 'expression', list( bquote( kappa == .(kappa_var)))))
# legend("topleft", legend = as.expression(sapply(kappa_var, function(var) bquote(bold("kappa") == .(var)))))
legend("topleft", legend = as.expression(sapply(kappa_var, function(var) bquote("kappa" == bold(.(var))))))
# legend("topleft", legend = as.expression(sapply(kappa_var, function(var) paste0(bold("kappa"), "==", var))))
dev.off()
