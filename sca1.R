#!/usr/bin/env Rscript

# interesting way to present.
x <- matrix(1:10, ncol = 2)

# So the scale() function will actually not scale if you don't want.
# which is funny, but that's because it does something else:
(centered.x <- scale(x, scale = FALSE))
# yes, it just centres in the above case, but that's interesting to look at.

# of course, it CAN scale ...
centered.scaled.x <- scale(x)

# now we'll use cov() to "show up" that after scaling
# both vectors are actually the same.
co <- cov(centered.scaled.x)
