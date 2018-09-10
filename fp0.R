#!/usr/bin/env Rscript
# experiments with R and floating point.

print("Infamously, floating points at high precisions can be a bit of a problem.")
print("And a pretty unavoidable problem at that due to underlying binary nature of")
print("computers.")

print("Let's take 0.2, and print it out with the (maximum) 22 dec pt precision")

c <- 0.2
print(c, digits=22)
# without quotes, this will be a numeric, which defautls to double.
# i.e. as.double() would have no effect on this.
# c <- "0.2000000000000000111022"
print("In fact withe sprintf() func, we can squeeze more digits out")
cs <- sprintf("%.28f", c)
print("However all this will be a bit pointless because alot of computations don't go that far.")
cc2 <- "0.200000000000000013"
cc2l <- nchar(cc2)
po <- sprintf("take this %i decptprecision characterstr: %s", cc2l-2, cc2)
po
po <- sprintf("as.double it and compare with 0.2: %i (1 for TRUE)", as.double(cc2) == c)
po
cc2 <- "0.20000000000000002"
cc2l <- nchar(cc2)
po <- sprintf("Now, ... take this %i decptprecision characterstr: %s", cc2l-2, cc2)
po
po <- sprintf("as.double it and compare with 0.2: %i (1 for TRUE)", as.double(cc2) == c)
po
po <- sprintf("So at what point does it become different?")
po
cc2 <- "0.20000000000000003"
cc2l <- nchar(cc2)
po <- sprintf("Now, ... take this %i decptprecision characterstr: %s", cc2l-2, cc2)
po
po <- sprintf("as.double it and compare with 0.2: %i (1 for TRUE)", as.double(cc2) == c)
po
print("Finally, we hve a different number by changing 17th precision digit")
print("but we had to increase by 2, 1 didn't manage it.")
