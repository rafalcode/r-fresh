#!/usr/bin/env Rscript
# ddm degenerate design matrices from:
# https://rstudio-pubs-static.s3.amazonaws.com/6311_a09169ad892f4f5499874751a5fa822d.html
# did a search on some bits of this, couldn't find who wrote it
# this is a short form of the main ddm.R
library(ggplot2)
library(Cairo)
library(Matrix)

# Degenerate design matrices
# One problem that arises sometimes in regression and related problems is a degenerate design matrix – or in other words perfect multicollinearity.

set.seed(101)
d <- data.frame(y = rnorm(5), x1 = 1:5, x2 = 1:5,
                f1 = rep(c("a", "b"), c(2, 3)),
                f2 = rep(c("A", "B"), c(3, 2)))

# Sometimes the problem is forehead-slappingly obvious, but what if it's not? You can construct your own design matrix using R's model.matrix() function, with the same formula (you only need the right-hand side) that you put into your model:

X1 <- model.matrix(~x1 + x2, data = d)
(r <- rankMatrix(X1))

## [1] 2
## attr(,"method")
## [1] "tolNorm2"
## attr(,"useGrad")
## [1] FALSE
## attr(,"tol")
## [1] 1.186e-14
# if r < ncol(X1)  ## this signals trouble
## [1] TRUE


# You can see which particular columns are involved by doing a singular value decomposition:

(s <- svd(X1))
# and seeing which columns are associated with the near-zero elements of s$d:
mcols <- s$v[, (s$d < .Machine$double.eps)]
setNames(mcols, colnames(X1))

## (Intercept)          x1          x2 
##      0.0000      0.7071     -0.7071

# This tells us that columns 2 and 3 of the design matrix (associated with x1 and x2 are perfectly collinear.
# Unfortunately I don't get this. These numbers are just the third column 

# A similar situation can arise with factors, especially when one tries to fit an interaction term where not all combinations of levels have been measured (either because they are logically impossible – structural zeros – or because they simply weren't measured.

# For example, the a:B combination is missing in this fake data set:

with(d, table(f1, f2))

##    f2
## f1  A B
##   a 2 0
##   b 1 2

# We will get similar results as when we tried to fit above:

lm10 <- lm(y ~ f1 * f2, data = d)

## 
## Call:
## lm(formula = y ~ f1 * f2, data = d)
## 
## Coefficients:
## (Intercept)          f1b          f2B      f1b:f2B  
##       0.113       -0.788        0.938           NA

library("nlme")
gls(y ~ f1 * f2, data = d)

## Error: computed "gls" fit is singular, rank 4

X2 <- model.matrix(~f1 * f2, data = d)
setNames(zapsmall(svd(X2)$v[, 4]), colnames(X2))

## (Intercept)         f1b         f2B     f1b:f2B 
##      0.0000      0.0000      0.7071     -0.7071

# As pointed out in the GLMM faq, the way to deal with this problem is to convert the problem into a one-way design (i.e. compute the interaction explicitly and drop unused levels):

d <- transform(d, f12 = droplevels(interaction(f1, f2)))
gls(y ~ f12, data = d)

## Generalized least squares fit by REML
##   Model: y ~ f12 
##   Data: d 
##   Log-restricted-likelihood: -1.898
## 
## Coefficients:
## (Intercept)      f12b.A      f12b.B 
##      0.1132     -0.7882      0.1494 
## 
## Degrees of freedom: 5 total; 2 residual
## Residual standard error: 0.4419

detach("package:nlme")

