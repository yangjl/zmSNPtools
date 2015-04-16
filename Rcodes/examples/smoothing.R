## how to do the smoothing curves!

#lo <- loess(gp[,2] ~ gp[,4])
#lines(predict(lo), col='red', lwd=2)
#You can set a smoothing parameter (typically between 0 and 1) here
smoothingSpline = smooth.spline(x, y, spar=1)
lines(smoothingSpline, col='red', lwd=2)

### fit linear line
reg1 <- lm(write~read)
par(cex=.8)
plot(read, write)
abline(reg1)




