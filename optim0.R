# https://stackoverflow.com/questions/15882323/r-robust-fitting-of-data-points-to-a-gaussian-function

# fitG =
# function(x,y,mu,sig,scale){
# 
#   f = function(p){
#     d = p[3]*dnorm(x,mean=p[1],sd=p[2])
#     sum((d-y)^2)
#   }
# 
#   optim(c(mu,sig,scale),f)
#  }

lossf <- function(p){
  d = p[3]*dnorm(x,mean=p[1],sd=p[2])
  sum((d-y)^2)
}


n <- 20
mu <- 50
sd <- 10
scale <- 1
da <- rnorm(n, mu, sd)

# optim(c(mu,sig,scale),f)

