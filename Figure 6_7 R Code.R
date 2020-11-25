setwd(dirname(rstudioapi::getActiveDocumentContext()$path))

library(stats)
library(rcosmo)
library(sp)
library(colf)
###################### Renyi and f_alpha functions #########################
#The function fRen computes the Renyi function for a chosen window of sky data(e.g.:- whole sky, large, medium windows,etc.).
fRen1 <-
  function(cmbdf,
           q.min = 1.01,
           q.max = 10,
           N = 20,
           k.box = log2(nside(cmbdf)) - 3,
           intensities = "I"){
    if (!is.CMBDataFrame(cmbdf)){
      stop("Argument must be a CMBDataFrame")}
    ns1 <- nside(cmbdf)
    pixind <- pix(cmbdf)
    nagrpix <- setdiff(1:(12 * ns1 ^ 2), pixind)
    field.comp <- rep(0, 12 * ns1 ^ 2)
    field.in <- cmbdf[, "I", drop = T]
    field.in <- field.in - minint
    field.final <- replace(field.comp, pixind, field.in)
    res.max <- log2(ns1)
    k.box <- log2(nside(cmbdf)) - 3
    npix <- 12 * 4 ^ k.box
    lev.diff <- 4 ^ (res.max - k.box)
    if (res.max - k.box > 0){
      nagrpix <- unique(ancestor(nagrpix, res.max - k.box))}
    agrpix <- setdiff((1:npix), nagrpix)
    mu <- vector(mode = "numeric", length = length(agrpix))
    field.total <- 0
    i <- 1
    for (j in agrpix){
      pixd <- (lev.diff * (j - 1) + 1):(lev.diff * j)
      field.total <- field.total + sum(field.final[pixd])
      mu[i] <- sum(field.final[pixd])
      i <- i + 1}
    mu <- mu / field.total
    Q <- seq(q.min, q.max, length.out = N)
    Tq <- vector(mode = "numeric", length = N)
    delta <- (1 / length(agrpix))
    ri <- 1
    for (q in Q){
      Tq[ri] <-  log2(sum(mu ^ q)) / log2(delta)
      ri <- ri + 1}
    Tqf <- data.frame(q = Q, tq = Tq)
    return(Tqf)}


#The function alp1 computes the alpha function for a chosen window of sky data(e.g.:- whole sky, large, medium windows,etc.).
alp1 <-
  function(cmbdf,
           q.min = 1.01,
           q.max = 10,
           N = 20,
           k.box = log2(nside(cmbdf)) - 3,
           intensities = "I"){
    if (!is.CMBDataFrame(cmbdf)){
      stop("Argument must be a CMBDataFrame")}
    ns1 <- nside(cmbdf)
    pixind <- pix(cmbdf)
    nagrpix <- setdiff(1:(12 * ns1 ^ 2), pixind)
    field.comp <- rep(0, 12 * ns1 ^ 2)
    field.in <- cmbdf[, "I", drop = T]
    field.in <- field.in - minint
    field.final <- replace(field.comp, pixind, field.in)
    res.max <- log2(ns1)
    k.box <- log2(nside(cmbdf)) - 3
    npix <- 12 * 4 ^ k.box
    lev.diff <- 4 ^ (res.max - k.box)
    if (res.max - k.box > 0){
      nagrpix <- unique(ancestor(nagrpix, res.max - k.box))}
    agrpix <- setdiff((1:npix), nagrpix)
    mu <- vector(mode = "numeric", length = length(agrpix))
    field.total <- 0
    i <- 1
    for (j in agrpix){
      pixd <- (lev.diff * (j - 1) + 1):(lev.diff * j)
      field.total <- field.total + sum(field.final[pixd])
      mu[i] <- sum(field.final[pixd])
      i <- i + 1}
    mu <- mu / field.total
    Q <- seq(q.min, q.max, length.out = N)
    alp <- vector(mode = "numeric", length = N)
    delta <- (1 / length(agrpix))
    ri <- 1
    for (q in Q){
      alp [ri] <- sum((mu ^ q / sum(mu ^ q)) * log2(mu)) / log2(delta)
      ri <- ri + 1}
    alp0  <- data.frame(q = Q, tq = alp)
    return(alp0)}


#The function falp1 computes the falpha function for a chosen window of sky data(e.g.:- whole sky, large, medium windows,etc.).
falp1 <-
  function(cmbdf,
           q.min = 1.01,
           q.max = 10,
           N = 20,
           k.box = log2(nside(cmbdf)) - 3,
           intensities = "I"){
    if (!is.CMBDataFrame(cmbdf)){
      stop("Argument must be a CMBDataFrame")}
    ns1 <- nside(cmbdf)
    pixind <- pix(cmbdf)
    nagrpix <- setdiff(1:(12 * ns1 ^ 2), pixind)
    field.comp <- rep(0, 12 * ns1 ^ 2)
    field.in <- cmbdf[, "I", drop = T]
    field.in <- field.in - minint
    field.final <- replace(field.comp, pixind, field.in)
    res.max <- log2(ns1)
    k.box <- log2(nside(cmbdf)) - 3
    npix <- 12 * 4 ^ k.box
    lev.diff <- 4 ^ (res.max - k.box)
    if (res.max - k.box > 0){
      nagrpix <- unique(ancestor(nagrpix, res.max - k.box))}
    agrpix <- setdiff((1:npix), nagrpix)
    mu <- vector(mode = "numeric", length = length(agrpix))
    field.total <- 0
    i <- 1
    for (j in agrpix){
      pixd <- (lev.diff * (j - 1) + 1):(lev.diff * j)
      field.total <- field.total + sum(field.final[pixd])
      mu[i] <- sum(field.final[pixd])
      i <- i + 1}
    mu <- mu / field.total
    Q <- seq(q.min, q.max, length.out = N)
    falp <- vector(mode = "numeric", length = N)
    delta <- (1 / length(agrpix))
    ri <- 1
    for (q in Q){
      falp [ri] <- sum((mu ^ q / sum(mu ^ q)) * log2(mu ^ q / sum(mu ^ q))) /
        log2(delta)
      ri <- ri + 1}
    falp0  <- data.frame(q = Q, tq = falp)
    return(falp0)}


############ Investigating Renyi functions #############
#Here, the behaviour of the Renyi function is being investigated to check whether the cosmic microwave background data are
#multifractal or not.

####################### For the whole sky ####################################

cmbdf <- CMBDataFrame("CMB_map_smica1024.fits")
minint <- min(cmbdf[, "I", drop = T])
q.min <- 1.01
q.max <- 2.0
N <- 20

########### Figure 6a-Sample Renyi function with the linear function
#Figure 6a gives the plot of the Sample Renyi function with the linear function for whole sky data.
Tq <- fRen1(cmbdf, q.min, q.max, N)
plot(
  Tq[, 1],
  Tq[, 2],
  ylab = "T(q)",
  xlab = "q",
  pch = 20,
  col = "blue",
  cex.main = 1.25,
  cex.lab = 1.25,
  cex.axis = 1)
segments(Tq[1, 1], Tq[1, 2], Tq[20, 1], Tq[20, 2], lwd = (2), col = "red")

Tq[, 3] <-
  (Tq[, 1] * (Tq[1, 2] - Tq[N, 2]) + Tq[N, 2] * Tq[1, 1] - Tq[1, 2] * Tq[N, 1]) /
  (Tq[1, 1] - Tq[N, 1])

########### Figure 6b-Difference between the Sample Renyi function and linear function
#Figure 6b gives the plot of the difference between the Sample Renyi function and linear function for whole sky data.
plot(
  Tq[, 1],
  Tq[, 2] - Tq[, 3],
  ylab = "Difference",
  xlab = "q",
  pch = 20,
  col = "blue")

x <- Tq[, 1]
b <- rep(1, 20)
y <- (Tq[, 2] - x + b)
Ren1 <- data.frame(x, y)

QM1 <- colf_nls(y ~ 0 + I(-(x ^ 2) + x), data = Ren1, lower = c(0))
#Coef(QM1) gives the estimated parameter a
coef(QM1)
Q11 <- (fitted(QM1) + x - 1)

########### Figure 6e-Sample Renyi function with the fitted log-normal model
#Figure 6e gives the plot of the Sample Renyi function with the fitted log-normal model for whole sky data.
plot(
  Tq[, 1],
  Tq[, 2],
  ylab = "T(q)",
  xlab = "q",
  pch = 20,
  col = "blue")
lines(Tq[, 1], Q11, col = "red", lwd = 2)

########### Figure 6f-Difference between the Sample Renyi function and the fitted log-normal model
#Figure 6f gives the plot of the difference between the Sample Renyi function and the fitted log-normal model for whole sky data.
plot(
  Tq[, 1],
  Q11 - Tq[, 2],
  ylab = "Difference",
  xlab = "q",
  pch = 20,
  col = "blue")

############# Root Mean Square Error (RMSE)
residuals <- Q11 - Tq[, 2]
sqrt(mean(residuals ^ 2))

q.min <- -10.0
q.max <- 10.0

#Computes the alpha function for whole sky data
Alp <- alp1(cmbdf, q.min, q.max, N)

########### Figure 6c-expression(paste("Sample ",  alpha, " versus ", q))
#Figure 6c gives the plot of the function alpha versus q for whole sky data.
plot(
  Alp[, 1],
  Alp[, 2],
  ylab = expression(paste(alpha(q))),
  xlab = "q",
  pch = 20,
  col = "red",
  cex.main = 1.25,
  cex.lab = 1.25,
  cex.axis = 1,
  type = "l")
points(
  Alp[, 1],
  Alp[, 2],
  ylab = expression(paste(alpha(q))),
  xlab = "q",
  pch = 19,
  col = "blue",
  cex.main = 1.25,
  cex.lab = 1.25,
  cex.axis = 1)

##### alpha interval

min(Alp[, 2])
max(Alp[, 2])
max(Alp[, 2]) - min(Alp[, 2])

#Computes the falpha function for whole sky data
Fq <- falp1(cmbdf, q.min, q.max, N)
plot(
  Fq[, 1],
  Fq[, 2],
  ylab = expression(paste(f[alpha](q))),
  xlab = "q",
  main = expression(paste("Sample ", f[alpha], " function")),
  pch = 20,
  col = "blue",
  cex.main = 1.25,
  cex.lab = 1.25,
  cex.axis = 1)

########### Figure 6d-expression(paste("Sample ", f(alpha), " vs ", alpha))
#Figure 6d gives the plot of the function falpha versus alpha for whole sky data.
plot(
  Alp[, 2],
  Fq[, 2],
  ylab = expression(paste(f(alpha))),
  xlab = expression(paste(alpha)),
  pch = 20,
  col = "red",
  cex.main = 1.25,
  cex.lab = 1.25,
  cex.axis = 1,
  type = "l")
points(
  Alp[, 2],
  Fq[, 2],
  ylab = expression(paste(f(alpha))),
  xlab = expression(paste(alpha)),
  pch = 19,
  col = "blue",
  cex.main = 1.25,
  cex.lab = 1.25,
  cex.axis = 1)

####################### For a large Window near the pole of the sphere ####################################

q.min <- 1.01
q.max <- 2.0
N <- 20
#Here, a large window of CMB sky data is chosen such that the area=1.231.
win <-
  CMBWindow(theta = c(3 * pi / 6, 3 * pi / 6, pi / 4, pi / 4),
            phi = c(0, pi / 2, pi / 2, 0))
cmbdf1 <- window(cmbdf, new.window = win)
Tq <- fRen1(cmbdf1, q.min, q.max, N)
#This figure gives the plot of the Sample Renyi function with the linear function for a large window of sky data.
plot(
  Tq[, 1],
  Tq[, 2],
  ylab = "T(q)",
  xlab = "q",
  pch = 20,
  col = "blue",
  cex.main = 1.25,
  cex.lab = 1.25,
  cex.axis = 1)
segments(Tq[1, 1], Tq[1, 2], Tq[20, 1], Tq[20, 2], lwd = (2), col = "red")

Tq[, 3] <-
  (Tq[, 1] * (Tq[1, 2] - Tq[N, 2]) + Tq[N, 2] * Tq[1, 1] - Tq[1, 2] * Tq[N, 1]) /
  (Tq[1, 1] - Tq[N, 1])

########### Figure 7b-Difference between the Sample Renyi function and linear function
#Figure 7b gives the plot of the difference between the Sample Renyi function and linear function for a large window of sky data.
plot(
  Tq[, 1],
  Tq[, 2] - Tq[, 3],
  ylab = "Difference",
  xlab = "q",
  pch = 20,
  col = "blue")

x <- Tq[, 1]
b <- rep(1, 20)
y <- (Tq[, 2] - x + b)
Ren1 <- data.frame(x, y)

QM1 <- colf_nls(y ~ 0 + I(-(x ^ 2) + x), data = Ren1, lower = c(0))
#Coef(QM1) gives the estimated parameter a
coef(QM1)
Q11 <- (fitted(QM1) + x - 1)

#This figure gives the plot of the Sample Renyi function with the fitted log-normal model for a large window of sky data.
plot(
  Tq[, 1],
  Tq[, 2],
  ylab = "T(q)",
  xlab = "q",
  pch = 20,
  col = "blue")
lines(Tq[, 1], Q11, col = "red", lwd = 2)

########### Figure 7c-Difference between the Sample Renyi function and the fitted log-normal model
#Figure 7c gives the plot of the difference between the Sample Renyi function and the fitted log-normal model for a large window of sky data.
plot(
  Tq[, 1],
  Q11 - Tq[, 2],
  ylab = "Difference",
  xlab = "q",
  pch = 20,
  col = "blue")

residuals <- Q11 - Tq[, 2]
sqrt(mean(residuals ^ 2))

q.min <- -10.0
q.max <- 10.0

#Computes the alpha function for a large window of sky data.
Alp <- alp1(cmbdf1, q.min, q.max, N)
#This figure gives the plot of the function alpha versus q for a large window of sky data.
plot(
  Alp[, 1],
  Alp[, 2],
  ylab = expression(paste(alpha(q))),
  xlab = "q",
  pch = 20,
  col = "blue",
  cex.main = 1.25,
  cex.lab = 1.25,
  cex.axis = 1)

min(Alp[, 2])
max(Alp[, 2])
max(Alp[, 2]) - min(Alp[, 2])

#Computes the falpha function for a large window of sky data.
Fq <- falp1(cmbdf1, q.min, q.max, N)
plot(
  Fq[, 1],
  Fq[, 2],
  ylab = expression(paste(f[alpha](q))),
  xlab = "q",
  pch = 20,
  col = "blue",
  cex.main = 1.25,
  cex.lab = 1.25,
  cex.axis = 1)

########### Figure 7a-expression(paste("Sample ", f(alpha), " vs ", alpha))
#Figure 7a gives the plot of the function falpha versus alpha for a large window of sky data.
plot(
  Alp[, 2],
  Fq[, 2],
  ylab = expression(paste(f(alpha))),
  xlab = expression(paste(alpha)),
  pch = 20,
  col = "red",
  cex.main = 1.25,
  cex.lab = 1.25,
  cex.axis = 1,
  type = "l")
points(
  Alp[, 2],
  Fq[, 2],
  ylab = expression(paste(f(alpha))),
  xlab = expression(paste(alpha)),
  pch = 19,
  col = "blue",
  cex.main = 1.25,
  cex.lab = 1.25,
  cex.axis = 1)

####################### For a small Window near the pole of the sphere ##########

q.min <- 1.01
q.max <- 2.0
N <- 20
#Here, a small window of CMB sky data is chosen such that the area=0.0596
win <-
  CMBWindow(theta = c(pi / 6, pi / 6, pi / 12, pi / 12),
            phi = c(0, pi / 5, pi / 5, 0))
cmbdf4 <- window(cmbdf, new.window = win)
Tq <- fRen1(cmbdf4, q.min, q.max, N)
#This figure gives the plot of the Sample Renyi function with the linear function for a small window of sky data.
plot(
  Tq[, 1],
  Tq[, 2],
  ylab = "T(q)",
  xlab = "q",
  pch = 20,
  col = "blue",
  cex.main = 1.25,
  cex.lab = 1.25,
  cex.axis = 1)
segments(Tq[1, 1], Tq[1, 2], Tq[20, 1], Tq[20, 2], lwd = (2), col = "red")

Tq[, 3] <-
  (Tq[, 1] * (Tq[1, 2] - Tq[N, 2]) + Tq[N, 2] * Tq[1, 1] - Tq[1, 2] * Tq[N, 1]) /
  (Tq[1, 1] - Tq[N, 1])

########### Figure 7e-Difference between the Sample Renyi function and linear function
#Figure 7e gives the plot of the difference between the Sample Renyi function and linear function for a small window of sky data.
plot(
  Tq[, 1],
  Tq[, 2] - Tq[, 3],
  ylab = "Difference",
  xlab = "q",
  pch = 20,
  col = "blue")

x <- Tq[, 1]
b <- rep(1, 20)
y <- (Tq[, 2] - x + b)
Ren1 <- data.frame(x, y)

QM1 <- colf_nls(y ~ 0 + I(-(x ^ 2) + x), data = Ren1, lower = c(0))
#Coef(QM1) gives the estimated parameter a
coef(QM1)
Q11 <- (fitted(QM1) + x - 1)

#This figure gives the plot of the Sample Renyi function with the fitted log-normal model for a small window of sky data.
plot(
  Tq[, 1],
  Tq[, 2],
  ylab = "T(q)",
  xlab = "q",
  pch = 20,
  col = "blue")
lines(Tq[, 1], Q11, col = "red", lwd = 2)

########### Figure 7f-Difference between the Sample Renyi function and the fitted log-normal model
#Figure 7f gives the plot of the difference between the Sample Renyi function and the fitted log-normal model for a small window of sky data.
plot(
  Tq[, 1],
  Q11 - Tq[, 2],
  ylab = "Difference",
  xlab = "q",
  pch = 20,
  col = "blue")

residuals <- Q11 - Tq[, 2]
sqrt(mean(residuals ^ 2))

q.min <- -10.0
q.max <- 10.0

#Computes the alpha function for a small window of sky data.
Alp <- alp1(cmbdf4, q.min, q.max, N)
#This figure gives the plot of the function alpha versus q for a small window of sky data.
plot(
  Alp[, 1],
  Alp[, 2],
  ylab = expression(paste(alpha(q))),
  xlab = "q",
  pch = 20,
  col = "blue",
  cex.main = 1.25,
  cex.lab = 1.25,
  cex.axis = 1)

min(Alp[, 2])
max(Alp[, 2])
max(Alp[, 2]) - min(Alp[, 2])

#Computes the falpha function for a small window of sky data.
Fq <- falp1(cmbdf4, q.min, q.max, N)
plot(
  Fq[, 1],
  Fq[, 2],
  ylab = expression(paste(f[alpha](q))),
  xlab = "q",
  pch = 20,
  col = "blue",
  cex.main = 1.25,
  cex.lab = 1.25,
  cex.axis = 1)

########### Figure 7d-expression(paste("Sample ", f(alpha), " vs ", alpha))
#Figure 7d gives the plot of the function falpha versus alpha for a small window of sky data.
plot(
  Alp[, 2],
  Fq[, 2],
  ylab = expression(paste(f(alpha))),
  xlab = expression(paste(alpha)),
  pch = 20,
  col = "red",
  cex.main = 1.25,
  cex.lab = 1.25,
  cex.axis = 1,
  type = "l")
points(
  Alp[, 2],
  Fq[, 2],
  ylab = expression(paste(f(alpha))),
  xlab = expression(paste(alpha)),
  pch = 19,
  col = "blue",
  cex.main = 1.25,
  cex.lab = 1.25,
  cex.axis = 1)

####################### For a medium Window near the pole of the sphere ####################################

q.min <- 1.01
q.max <- 2.0
N <- 20
#Here, a medium window of CMB sky data is chosen such that the area=0.4056 
win <-
  CMBWindow(theta = c(pi / 3.5, pi / 3.5, pi / 10, pi / 10),
            phi = c(0, pi / 2, pi / 2, 0))
cmbdf2 <- window(cmbdf, new.window = win)
Tq <- fRen1(cmbdf2, q.min, q.max, N)
#This figure gives the plot of the Sample Renyi function with the linear function for a medium window of sky data.
plot(
  Tq[, 1],
  Tq[, 2],
  ylab = "T(q)",
  xlab = "q",
  pch = 20,
  col = "blue",
  cex.main = 1.25,
  cex.lab = 1.25,
  cex.axis = 1)
segments(Tq[1, 1], Tq[1, 2], Tq[20, 1], Tq[20, 2], lwd = (2), col = "red")

Tq[, 3] <-
  (Tq[, 1] * (Tq[1, 2] - Tq[N, 2]) + Tq[N, 2] * Tq[1, 1] - Tq[1, 2] * Tq[N, 1]) /
  (Tq[1, 1] - Tq[N, 1])

#This figure gives the plot of the difference between the Sample Renyi function and linear function for a medium window of sky data.
plot(
  Tq[, 1],
  Tq[, 2] - Tq[, 3],
  ylab = "difference",
  xlab = "q",
  pch = 20,
  col = "blue")

x <- Tq[, 1]
b <- rep(1, 20)
y <- (Tq[, 2] - x + b)
Ren1 <- data.frame(x, y)

QM1 <- colf_nls(y ~ 0 + I(-(x ^ 2) + x), data = Ren1, lower = c(0))
#Coef(QM1) gives the estimated parameter a
coef(QM1)
Q11 <- (fitted(QM1) + x - 1)

#This figure gives the plot of the Sample Renyi function with the fitted log-normal model for a medium window of sky data.
plot(
  Tq[, 1],
  Tq[, 2],
  ylab = "T(q)",
  xlab = "q",
  pch = 20,
  col = "blue")
lines(Tq[, 1], Q11, col = "red", lwd = 2)

residuals <- Q11 - Tq[, 2]
sqrt(mean(residuals ^ 2))

q.min <- -10.0
q.max <- 10.0

#Computes the alpha function for a medium window of sky data.
Alp <- alp1(cmbdf2, q.min, q.max, N)
#This figure gives the plot of the function alpha versus q for a medium window of sky data.
plot(
  Alp[, 1],
  Alp[, 2],
  ylab = expression(paste(alpha(q))),
  xlab = "q",
  pch = 20,
  col = "blue",
  cex.main = 1.25,
  cex.lab = 1.25,
  cex.axis = 1)

min(Alp[, 2])
max(Alp[, 2])
max(Alp[, 2]) - min(Alp[, 2])

#Computes the falpha function for a medium window of sky data.
Fq <- falp1(cmbdf2, q.min, q.max, N)
plot(
  Fq[, 1],
  Fq[, 2],
  ylab = expression(paste(f[alpha](q))),
  xlab = "q",
  pch = 20,
  col = "blue",
  cex.main = 1.25,
  cex.lab = 1.25,
  cex.axis = 1)

#This figure gives the plot of the function falpha versus alpha for a medium window of sky data.
plot(
  Alp[, 2],
  Fq[, 2],
  ylab = expression(paste(f(alpha))),
  xlab = expression(paste(alpha)),
  pch = 20,
  col = "red",
  cex.main = 1.25,
  cex.lab = 1.25,
  cex.axis = 1,
  type = "l")
points(
  Alp[, 2],
  Fq[, 2],
  ylab = expression(paste(f(alpha))),
  xlab = expression(paste(alpha)),
  pch = 19,
  col = "blue",
  cex.main = 1.25,
  cex.lab = 1.25,
  cex.axis = 1)

###################### For a very small window near the pole of the sphere #########################

q.min <- 1.01
q.max <- 2.0
N <- 20
#Here, a very small window of CMB sky data is chosen such that the area=0.0017
win <-
  CMBWindow(theta = c(pi / 15, pi / 15, pi / 20, pi / 20),
            phi = c(0, pi / 18, pi / 18, 0))
cmbdf9 <- window(cmbdf, new.window = win)
Tq <- fRen1(cmbdf9, q.min, q.max, N)
#This figure gives the plot of the Sample Renyi function with the linear function for a very small window of sky data.
plot(
  Tq[, 1],
  Tq[, 2],
  ylab = "T(q)",
  xlab = "q",
  pch = 20,
  col = "blue",
  cex.main = 1.25,
  cex.lab = 1.25,
  cex.axis = 1)
segments(Tq[1, 1], Tq[1, 2], Tq[20, 1], Tq[20, 2], lwd = (2), col = "red")

Tq[, 3] <-
  (Tq[, 1] * (Tq[1, 2] - Tq[N, 2]) + Tq[N, 2] * Tq[1, 1] - Tq[1, 2] * Tq[N, 1]) /
  (Tq[1, 1] - Tq[N, 1])

#This figure gives the plot of the difference between the Sample Renyi function and linear function for a very small window of sky data.
plot(
  Tq[, 1],
  Tq[, 2] - Tq[, 3],
  ylab = "difference",
  xlab = "q",
  pch = 20,
  col = "blue")

x <- Tq[, 1]
b <- rep(1, 20)
y <- (Tq[, 2] - x + b)
Ren1 <- data.frame(x, y)

QM1 <- colf_nls(y ~ 0 + I(-(x ^ 2) + x), data = Ren1, lower = c(0))
#Coef(QM1) gives the estimated parameter a
coef(QM1)
Q11 <- (fitted(QM1) + x - 1)

#This figure gives the plot of the Sample Renyi function with the fitted log-normal model for a very small window of sky data.
plot(
  Tq[, 1],
  Tq[, 2],
  ylab = "T(q)",
  xlab = "q",
  pch = 20,
  col = "blue")
lines(Tq[, 1], Q11, col = "red", lwd = 2)

residuals <- Q11 - Tq[, 2]
sqrt(mean(residuals ^ 2))

q.min <- -10.0
q.max <- 10.0

#Computes the alpha function for a very small window of sky data.
Alp <- alp1(cmbdf9, q.min, q.max, N)

#This figure gives the plot of the function alpha versus q for a very small window of sky data.
plot(
  Alp[, 1],
  Alp[, 2],
  ylab = expression(paste(alpha(q))),
  xlab = "q",
  pch = 20,
  col = "blue",
  cex.main = 1.25,
  cex.lab = 1.25,
  cex.axis = 1)

min(Alp[, 2])
max(Alp[, 2])
max(Alp[, 2]) - min(Alp[, 2])

#Computes the falpha function for a very small window of sky data.
Fq <- falp1(cmbdf9, q.min, q.max, N)
plot(
  Fq[, 1],
  Fq[, 2],
  ylab = expression(paste(f[alpha](q))),
  xlab = "q",
  pch = 20,
  col = "blue",
  cex.main = 1.25,
  cex.lab = 1.25,
  cex.axis = 1)

#This figure gives the plot of the function falpha versus alpha for a very small window of sky data.
plot(
  Alp[, 2],
  Fq[, 2],
  ylab = expression(paste(f(alpha))),
  xlab = expression(paste(alpha)),
  pch = 20,
  col = "red",
  cex.main = 1.25,
  cex.lab = 1.25,
  cex.axis = 1,
  type = "l")
points(
  Alp[, 2],
  Fq[, 2],
  ylab = expression(paste(f(alpha))),
  xlab = expression(paste(alpha)),
  pch = 19,
  col = "blue",
  cex.main = 1.25,
  cex.lab = 1.25,
  cex.axis = 1)

