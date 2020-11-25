library(sp)
library(RandomFieldsUtils)
library(RandomFields)
library(stats)
library(rcosmo)
library(colf)

###################### Renyi and f_alpha functions #########################
#The function fRen computes the Renyi function for a chosen window of sky data(e.g.:- whole sky, large, medium windows,etc.).

fRen1 <- function(cmbdf, q.min = 1.01, q.max = 10, N = 20, k.box = log2(nside(cmbdf)) - 3, intensities = "I")
{
  if (!is.CMBDataFrame(cmbdf)) {
    stop("Argument must be a CMBDataFrame")
  }
  ns1 <- nside(cmbdf)
  pixind <- pix(cmbdf)
  nagrpix <- setdiff(1:(12 * ns1^2), pixind)
  field.comp <- rep(0, 12 * ns1^2)
  field.in <- cmbdf[, "I", drop = T]
  field.in <- field.in - minint
  field.final <- replace(field.comp, pixind, field.in)
  res.max <- log2(ns1)
  k.box <- log2(nside(cmbdf)) - 3
  npix <- 12 * 4^k.box
  lev.diff <- 4^(res.max - k.box)
  if (res.max - k.box > 0) {
    nagrpix <- unique(ancestor(nagrpix, res.max - k.box))
  }
  agrpix <- setdiff((1:npix), nagrpix)
  mu <- vector(mode = "numeric", length = length(agrpix))
  field.total <- 0
  i <- 1
  for (j in agrpix) {
    pixd <- (lev.diff * (j - 1) + 1):(lev.diff * j)
    field.total <- field.total + sum(field.final[pixd])
    mu[i] <- sum(field.final[pixd])
    i <- i + 1
  }
  mu <- mu/field.total
  Q <- seq(q.min, q.max, length.out = N)
  Tq <- vector(mode = "numeric", length = N)
  delta <- (1/length(agrpix))
  ri <- 1
  for (q in Q) {
    Tq[ri] <-  log2(sum(mu^q))/log2(delta)
    ri <- ri + 1
  }
  Tqf <- data.frame(q = Q, tq = Tq)
  return(Tqf)
}


#The function alp1 computes the alpha function for a chosen window of sky data(e.g.:- whole sky, large, medium windows,etc.).
alp1 <- function (cmbdf, q.min = 1.01, q.max = 10, N = 20, k.box = log2(nside(cmbdf)) - 3, intensities = "I") 
{
  if (!is.CMBDataFrame(cmbdf)) {
    stop("Argument must be a CMBDataFrame")
  }
  ns1 <- nside(cmbdf)
  pixind <- pix(cmbdf)
  nagrpix <- setdiff(1:(12 * ns1^2), pixind)
  field.comp <- rep(0, 12 * ns1^2)
  field.in <- cmbdf[, "I", drop = T]
  field.in <- field.in - minint
  field.final <- replace(field.comp, pixind, field.in)
  res.max <- log2(ns1)
  k.box <- log2(nside(cmbdf)) - 3
  npix <- 12 * 4^k.box
  lev.diff <- 4^(res.max - k.box)
  if (res.max - k.box > 0) {
    nagrpix <- unique(ancestor(nagrpix, res.max - k.box))
  }
  agrpix <- setdiff((1:npix), nagrpix)
  mu <- vector(mode = "numeric", length = length(agrpix))
  field.total <- 0
  i <- 1
  for (j in agrpix) {
    pixd <- (lev.diff * (j - 1) + 1):(lev.diff * j)
    field.total <- field.total + sum(field.final[pixd])
    mu[i] <- sum(field.final[pixd])
    i <- i + 1
  }
  mu <- mu/field.total
  Q <- seq(q.min, q.max, length.out = N)
  alp <- vector(mode = "numeric", length = N)
  delta <- (1/length(agrpix))
  ri <- 1
  for (q in Q) {
    alp [ri] <- sum((mu^q/sum(mu^q))*log2(mu))/log2(delta)
    ri <- ri + 1
  }
  alp0  <- data.frame(q = Q, tq = alp )
  return(alp0)
}


#The function falp1 computes the falpha function for a chosen window of sky data(e.g.:- whole sky, large, medium windows,etc.).

falp1 <- function (cmbdf, q.min = 1.01, q.max = 10, N = 20, k.box = log2(nside(cmbdf)) - 3, intensities = "I") 
{
  if (!is.CMBDataFrame(cmbdf)) {
    stop("Argument must be a CMBDataFrame")
  }
  ns1 <- nside(cmbdf)
  pixind <- pix(cmbdf)
  nagrpix <- setdiff(1:(12 * ns1^2), pixind)
  field.comp <- rep(0, 12 * ns1^2)
  field.in <- cmbdf[, "I", drop = T]
  field.in <- field.in - minint
  field.final <- replace(field.comp, pixind, field.in)
  res.max <- log2(ns1)
  k.box <- log2(nside(cmbdf)) - 3
  npix <- 12 * 4^k.box
  lev.diff <- 4^(res.max - k.box)
  if (res.max - k.box > 0) {
    nagrpix <- unique(ancestor(nagrpix, res.max - k.box))
  }
  agrpix <- setdiff((1:npix), nagrpix)
  mu <- vector(mode = "numeric", length = length(agrpix))
  field.total <- 0
  i <- 1
  for (j in agrpix) {
    pixd <- (lev.diff * (j - 1) + 1):(lev.diff * j)
    field.total <- field.total + sum(field.final[pixd])
    mu[i] <- sum(field.final[pixd])
    i <- i + 1
  }
  mu <- mu/field.total
  Q <- seq(q.min, q.max, length.out = N)
  falp <- vector(mode = "numeric", length = N)
  delta <- (1/length(agrpix))
  ri <- 1
  for (q in Q) {
    falp [ri] <- sum((mu^q/sum(mu^q))*log2(mu^q/sum(mu^q)))/log2(delta)
    ri <- ri + 1
  }
  falp0  <- data.frame(q = Q, tq = falp )
  return(falp0)
}


############################ RMex
####################### Large Window ####################################
#Here, a multifractal random field is generated on a large window of sky data from a Gaussian mother field Y(x) with exponential
#covariance model and variance=2.

cmbdf <- CMBDataFrame(nside = 1024, I = rep(0, 12 * 1024 ^ 2))

q.min <- 0.5
q.max <- 3.0
N <- 40
win <-
  CMBWindow(theta = c(3 * pi / 6, 3 * pi / 6, pi / 4, pi / 4),
            phi = c(0, pi / 2, pi / 2, 0))
cmbdf2<- window(cmbdf, new.window = win)

df2 <- coords(cmbdf2, new.coords = "cartesian")

K <- 40
for (i in 0:K){
  model1 <- RMexp(var=2, scale=(3^i))
  f1 <- RFsimulate(x=df2$x, y=df2$y, z=df2$z, model=model1, spConform=FALSE)
  cmbdf2$I <- cmbdf2$I + f1
}
save.image(file = "SimMultField1.RData")
load("SimMultField1.RData")
cmbdf2$I <- exp(cmbdf2$I-K-1)

######### Figure 4(a)-This figure plots the realization of a multifractal random field generated from a Gaussian mother field Y(x) with 
#exponential covariance model and variance=2.
#The field is rescaled to match rcosmo colour pallet
cmbdf1 <- cmbdf2
cmbdf1$I <- (cmbdf2$I -mean(cmbdf2$I))*10^{15.6}
plot(cmbdf1,back.col="white",ylab="",xlab="",zlab="")

minint <- min(cmbdf2$I)
Tq <- fRen1(cmbdf2, q.min, q.max, N)
#This figure gives the plot of the Sample Renyi function with the linear function for the simulated large window of sky data.
plot(Tq[,1], Tq[,2], ylab ="T(q)", xlab = "q", main = "Sample Renyi function and linear function", pch = 20, col = "blue", cex.main=1.25, cex.lab=1.25, cex.axis=1)
segments(Tq[1,1], Tq[1,2], Tq[N,1], Tq[N,2], lwd = (2), col= "red")

Tq[,3] <- (Tq[,1]*(Tq[1,2]-Tq[N,2])+Tq[N,2]*Tq[1,1]-Tq[1,2]*Tq[N,1])/(Tq[1,1]-Tq[N,1])

#This figure gives the plot of the difference between the Sample Renyi function and linear function for the simulated large window of sky data.
plot(Tq[,1],Tq[,2]-Tq[,3], ylab ="difference", xlab = "q", main = "Difference of sample Renyi function and linear function", pch = 20, col = "blue")

x<-Tq[,1]
b<-rep(1,20)
y<-(Tq[,2]-x+b)
Ren1<-data.frame(x,y)

QM1<-colf_nls(y ~ 0+I(-(x^2)+x), data = Ren1, lower = c(0))

#Coef(QM1) gives the estimated parameter a
coef(QM1)
Q11<-(fitted(QM1)+x-1)

######### Figure 4(b)-This figure gives the plot of the Sample Renyi function with the fitted log-normal model for the simulated
#large window of sky data.

plot(Tq[,1],Tq[,2], ylab ="T(q)", xlab = "q", pch = 20, col = "blue")
lines(Tq[,1],Q11, col="red", lwd=2)

#This figure gives the plot of the difference between the Sample Renyi function and the fitted log-normal model for the simulated large window of sky data.
plot(Tq[,1],Q11-Tq[,2], ylab ="difference", xlab = "q", main = "Difference of Sample Renyi function and the fitted Log-Normal model", pch = 20, col = "blue")

residuals<- Q11-Tq[,2]
sqrt(mean(residuals^2))

q.min <- -10.0
q.max <- 10.0

Alp <- alp1(cmbdf2, q.min, q.max, N)

#This figure gives the plot of the function alpha versus q for the simulated large window of sky data.
plot(Alp[,1], Alp[,2], ylab =expression(paste(alpha(q))), xlab = "q", main = expression(paste("Sample ",  alpha, " versus ", q)), pch = 20, col = "blue", cex.main=1.25, cex.lab=1.25, cex.axis=1)

min(Alp[,2])
max(Alp[,2])
max(Alp[,2])- min(Alp[,2])

Fq <- falp1(cmbdf2, q.min, q.max, N)
plot(Fq[,1], Fq[,2], ylab =expression(paste(f[alpha](q))), xlab = "q", main = expression(paste("Sample ", f[alpha], " function")), pch = 20, col = "blue", cex.main=1.25, cex.lab=1.25, cex.axis=1)

########### Figure 4c-This figure gives the plot of the function falpha versus alpha for the simulated large window of sky data.
plot(Alp[,2], Fq[,2], ylab =expression(paste(f(alpha))), xlab = expression(paste(alpha)), pch = 20, col = "red", cex.main=1.25, cex.lab=1.25, cex.axis=1,type="l")
points(Alp[,2], Fq[,2], ylab =expression(paste(f(alpha))), xlab = expression(paste(alpha)), pch = 19, col = "blue", cex.main=1.25, cex.lab=1.25, cex.axis=1)




