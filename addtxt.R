#!/usr/bin/env Rscript
# using the text() commnd to add text on an R plot
library(Cairo)

K <- 16
CairoPNG("tx0.png", 800, 800)
plot(-1:1, -1:1, type = "n", xlab = "Re", ylab = "Im")
text(exp(1i * 2 * pi * (1:K) / K), col = 2)
dev.off()

## The following two examples use latin1 characters: these may not
## appear correctly (or be omitted entirely).
CairoPNG("tx1.png", 800, 800)
plot(1:20, 1:20, main = "text(...) examples\n~~~~~~~~~~~~~~", sub = "R is GNU ©, but not ® ...")
mtext("«Latin-1 accented chars»: éè øØ å<Å æ<Æ", side = 3)
points(c(6,2), c(2,1), pch = 3, cex = 4, col = "red")
text(6, 2, "the text is CENTERED around (x,y) = (6,2) by default", cex = .8)
text(2, 1, "or Left/Bottom - JUSTIFIED at (2,1) by 'adj = c(0,0)'", adj = c(0,0))
text(4, 9, expression(hat(beta) == (X^t * X)^{-1} * X^t * y))
text(4, 8.4, "expression(hat(beta) == (X^t * X)^{-1} * X^t * y)", cex = .75)
text(4, 7, expression(bar(x) == sum(frac(x[i], n), i==1, n)))

## Two more latin1 examples
text(5, 10.2, "Le français, c'est façile: Règles, Liberté, Egalité, Fraternité...")
text(5, 9.8, "Jetz no chli züritüütsch: (noch ein bißchen Zürcher deutsch)")
dev.off()
