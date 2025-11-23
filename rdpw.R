# How to read in classic pairwise distance data
# A simple format is space spearated labels on one line
# and on the second line the distances, also space separated.

# We'll also used that most simple of input functions readLines()
# so here is such a pairwise dataset of road distancees in miles from some cities in England and Scotland.
rl <- readLines("distsco.txt")

# So it's two lines and though this function is a ltitle raw
# you know what you're getting
# each line is a single member and with tw lines
# we get a 2 element vector because there's two lines.

# Now these two need to be parsed, but we know the numer of elements
# in each line will be different, so it does not lend itself to a vector of vectors
# At the same time, the vectros are linked in our analysis and there is
# a relationship between the two
# so it would be good to keep them together
# Good news! strsplit will automatically do this!
# as it will return a list (which is happy to have differently sized elements)
# of two vectors, each vector is accessed via the double square brackets method
# i.e. [[1]] and [[2]].
pwc <- strsplit(rl, " ") # apirwaise container.

# the realationship between the two szie should be as follows
# but beware the division by two could create a float
sz1 <- length(pwc[[1]])
nn <- sz1*(sz1-1)/2
sz2 <- length(pwc[[2]])

cat(paste0("     "))
nams <- sprintf("%5s", pwc[[1]][sz1:2])
cat(paste0(nams, collapse=" "))
cat(paste0("\n"))
cs <- cumsum((sz1-1):1)
j <- 1
k <- 1
for(i in cs) {
    # cat(paste0(j:i, collapse=" "))
    cat(paste0(pwc[[1]][k], " "))
    # cat(paste0(pwc[[2]][i:j], collapse=" "))
    vals <- sprintf("%5i", as.integer(pwc[[2]][i:j]))
    cat(paste0(vals, collapse=" "))
    cat(paste0("\n"))
    j <- i+1
    k <- k+1
}
