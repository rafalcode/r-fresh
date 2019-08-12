#!/usr/bin/Rscript
library(Rcpp)

cppsrcstr <-    'int g(int n)
                {
                    if (n < 2)
                        return(n);
                    return(g(n-1) + g(n-2));
                }'

print cppsrcstr
# cppFunction(cppsrcstr)
# sapply(0:10, g)
