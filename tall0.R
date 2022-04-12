# an invesitgation into tall format, sort of concatenate columns into one columns, one after the other.
# required by lattice arangement of figures.
library(lattice) # xyplot

# First generate 10 rows of time series data (daily timpe points)
t <- seq(as.Date('2009-01-01'),by='days',length=10)
X <- rnorm(10,0,1) # stdev 1, 2, and 4
Y <- rnorm(10,0,2)
Z <- rnorm(10,0,4)
dat <- data.frame(t,X,Y,Z)
# r$> ls()
# [1] "dat" "t"   "X"   "Y"   "Z"
# 
# r$> dat
#             t          X          Y          Z
# 1  2009-01-01  0.2116050 -0.2602909 -2.8374149
# 2  2009-01-02 -0.2373027  0.9238210 -1.5875828
# 3  2009-01-03 -1.0370440  2.8673773  2.1197380
# 4  2009-01-04 -1.3991553  0.7803139 -5.3262568
# 5  2009-01-05 -0.6413231 -5.2925434  3.2690545
# 6  2009-01-06  1.1815153  2.2780925  2.5784696
# 7  2009-01-07 -0.1917103 -4.7021496  1.5889211
# 8  2009-01-08  0.1929830 -0.6536213  2.4349112
# 9  2009-01-09  0.8692131  1.2498175 -0.2119724
# 10 2009-01-10 -1.0610302  2.7489315  1.5454169

# so it's the reshape() that will take care of this.
# seems to be a bit based on time series, with the "time" option appearing.
m <- reshape(dat, direction = "long", varying = 2:4, v.names = "price", idvar = "t", timevar = "symbol", times = names(dat)[2:4],   new.row.names = 1:30)
#             t symbol      price
# 1  2009-01-01      X  0.2116050
# 2  2009-01-02      X -0.2373027
# 3  2009-01-03      X -1.0370440
# 4  2009-01-04      X -1.3991553
# 5  2009-01-05      X -0.6413231
# 6  2009-01-06      X  1.1815153
# 7  2009-01-07      X -0.1917103
# 8  2009-01-08      X  0.1929830
# 9  2009-01-09      X  0.8692131
# 10 2009-01-10      X -1.0610302
# 11 2009-01-01      Y -0.2602909
# 12 2009-01-02      Y  0.9238210
# 13 2009-01-03      Y  2.8673773
# 14 2009-01-04      Y  0.7803139
# 15 2009-01-05      Y -5.2925434
# 16 2009-01-06      Y  2.2780925
# 17 2009-01-07      Y -4.7021496
# 18 2009-01-08      Y -0.6536213
# 19 2009-01-09      Y  1.2498175
# 20 2009-01-10      Y  2.7489315
# 21 2009-01-01      Z -2.8374149
# 22 2009-01-02      Z -1.5875828
# 23 2009-01-03      Z  2.1197380
# 24 2009-01-04      Z -5.3262568
# 25 2009-01-05      Z  3.2690545
# 26 2009-01-06      Z  2.5784696
# 27 2009-01-07      Z  1.5889211
# 28 2009-01-08      Z  2.4349112
# 29 2009-01-09      Z -0.2119724
# 30 2009-01-10      Z  1.5454169

png("tall0.png")
xyplot(price ~ t | symbol, data=m ,type ="l", layout = c(1,3) ) # a lattice func.
dev.off()
