#!/usr/bin/env Rscript
# adventures in smooth spline.
library(Cairo)

# the famous cars dataset, apparently 1920 is the year.
# note that distance is the stopping distance, so it's more a test of the brakes that the motor.
cars.sp0 <- with(cars, smooth.spline(dist))
cars.sp1 <- with(cars, smooth.spline(speed))
cars.spl <- with(cars, smooth.spline(speed, dist))

# Cairo image templatG
CairoPNG("smoo0.png", 1600, 800)
par(mfrow=c(1,2))
plot(cars$speed, main = "cars$speed")
plot(cars$dist, main = "cars$dist")
dev.off()

# Cairo image templatG
CairoPNG("smoo0.png", 1600, 800)
par(mfrow=c(1,2))
plot(cars$speed, main = "cars$speed & spline")
lines(cars.sp0, col = "skyblue")
plot(cars$dist, main = "cars$dist & spline")
lines(cars.sp1, col = "darkorange")
dev.off()

CairoPNG("smoo2.png", 800, 800)
plot(dist ~ speed, data = cars, main = "data(cars)  &  smoothing splines")
lines(cars.spl, col = "blue")
# put plot command here
dev.off()
