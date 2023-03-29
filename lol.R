#!/usr/bin/env Rscript
# this script does what?
# from:
# https://stackoverflow.com/questions/14848172/appending-a-list-to-a-list-of-lists-in-r
# Initial list:
lol <- list()
times <- 5
rnels <- 1+as.integer(runif(times)*3)

# Now the new experiments
for(i in 1:times){
  # lol[[length(lol)+1]] <- list(sample(1:3)) # example a bit inane.
  lol[[length(lol)+1]] <- sample(letters[22:26], rnels[i], replace=T)
}
write.csv(t(lol), "lol_.csv")
