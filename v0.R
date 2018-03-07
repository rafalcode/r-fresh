#!/usr/bin/Rscript
# #!/usr/bin/env Rscript
library(inline)
library(Rcpp)

src <- '
     Rcpp::NumericVector vec(vx);
     double p = Rcpp::as<double>(dd);
     double sum = 0.0;
     for (int i=0; i<vec.size(); i++) {
         sum += pow(vec[i], p);
     }
     return Rcpp::wrap(sum);'

fun <- cxxfunction(signature(vx="numeric", dd="numeric"), src, plugin="Rcpp")
fun(1:4,2)
fun(1:4,2.2)
