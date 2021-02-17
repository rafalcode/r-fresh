#!/usr/bin/env Rscript
# trying out bash expansions from R.

# note that brace expansion won;t really work as it's not really
# bash that the system call gets, closer to sh.
# truth is you might need bash -c, but this is complicated for a loop
# subshell will work though
# cmd <- paste0("for i in `seq 1 3`; do echo \"$i\";done")
# cmd <- paste0("for i in \"{1..4\"}; do echo \"my boy $i\"; done")
cmd <- paste0("for i in `seq 1 4`; do echo \"my boy $i\"; done")
cat(c(cmd, " si now executing.\n"))

outp <- system(cmd, intern=F)
cat(paste0("system coomand return value was ", outp, "\n"))

