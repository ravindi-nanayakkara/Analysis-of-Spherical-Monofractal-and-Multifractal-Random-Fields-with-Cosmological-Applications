setwd(dirname(rstudioapi::getActiveDocumentContext()$path))

library(stats)
library(rcosmo)
library(sp)

############## Renyi function
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
    for (j in agrpix) {
      pixd <- (lev.diff * (j - 1) + 1):(lev.diff * j)
      field.total <- field.total + sum(field.final[pixd])
      mu[i] <- sum(field.final[pixd])
      i <- i + 1}
    mu <- mu / field.total
    Q <- seq(q.min, q.max, length.out = N)
    Tq <- vector(mode = "numeric", length = N)
    delta <- (1 / length(agrpix))
    ri <- 1
    for (q in Q) {
      Tq[ri] <-  log2(sum(mu ^ q)) / log2(delta)
      ri <- ri + 1}
    Tqf <- data.frame(q = Q, tq = Tq)
    return(Tqf)}

#########################################################################

cmbdf <- CMBDataFrame("CMB_map_smica1024.fits")
minint <- min(cmbdf[, "I", drop = T])
q.min <- 1.01
q.max <- 2.0
N <- 20

####################### Model 4: Large Window near the pole of the sphere ####################################
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

x <- Tq[, 1]
b <- rep(1, 20)
y <- (Tq[, 2] - x + b)
Ren1 <- data.frame(x, y)

QM4 <-
  nls(
    y ~ A * (x + log2(gamma(x + (
      1 / 2
    ))) - 0.5 * log2(pi)),
    start = list(A = 0.2),
    data = Ren1,
    control = nls.control(
      maxiter = 1000,
      tol = 1e-02,
      minFactor = 1 / 1024,
      printEval = FALSE,
      warnOnly = TRUE),
    trace = TRUE)
#Coef(QM4) gives the estimated parameter A
coef(QM4)
Q14 <- (fitted(QM4) + x - 1)

#This figure gives the plot of the Sample Renyi function with the fitted Model 4 for a large window of sky data.
plot(
  Tq[, 1],
  Tq[, 2],
  ylab = "T(q)",
  xlab = "q",
  pch = 20,
  col = "blue")
lines(Tq[, 1], Q14, col = "red", lwd = 2)

########### Figure 8c-Difference between the Sample Renyi function and the fitted Model 4
#Figure 8c gives the plot of the difference between the Sample Renyi function and the fitted Model 4 for a large window of sky data.
plot(
  Tq[, 1],
  Q14 - Tq[, 2],
  ylab = "Difference",
  xlab = "q",
  pch = 20,
  col = "blue")

residuals <- Q14 - Tq[, 2]
sqrt(mean(residuals ^ 2))

########################## Model 5: Large Window near the pole of the sphere ##################################

QM5 <-
  nls(
    y ~ A * (B * x + log2(gamma(B * x + (
      1 / 2
    ))) - 0.5 * log2(pi)),
    start = list(A = 0.2, B = 2),
    data = Ren1,
    control = nls.control(
      maxiter = 10000000,
      tol = 1e-02,
      minFactor = 1 / 1024,
      printEval = FALSE,
      warnOnly = TRUE),
    trace = TRUE)
#Coef(QM5) gives the estimated parameters A and B.
coef(QM5)
Q15 <- (fitted(QM5) + x - 1)

#This figure gives the plot of the Sample Renyi function with the fitted Model 5 for a large window of sky data.
plot(
  Tq[, 1],
  Tq[, 2],
  ylab = "T(q)",
  xlab = "q",
  pch = 20,
  col = "blue")
lines(Tq[, 1], Q15, col = "red", lwd = 2)

########### Figure 8d-Difference between the Sample Renyi function and the fitted Model 5
#Figure 8d gives the plot of the difference between the Sample Renyi function and the fitted Model 5 for a large window of sky data.
plot(
  Tq[, 1],
  Q15 - Tq[, 2],
  ylab = "Difference",
  xlab = "q",
  pch = 20,
  col = "blue")

residuals <- Q15 - Tq[, 2]
sqrt(mean(residuals ^ 2))

####################### Model 6: Large Window near the pole of the sphere ####################################

QM6 <-
  nls(
    y ~ A * (x * log2(B) - x - log2(gamma(x + B)) + log2(gamma(B))),
    start = list(A = 0.2, B = 3),
    data = Ren1,
    control = nls.control(
      maxiter = 10000000,
      tol = 1e-02,
      minFactor = 1 / 1024,
      printEval = FALSE,
      warnOnly = TRUE),
    trace = TRUE,
    algorithm = "port",
    lower = c(0, 1),
    upper = c(100, 300))
#Coef(QM6) gives the estimated parameters A and B.
coef(QM6)
Q16 <- (fitted(QM6) + x - 1)

#This figure gives the plot of the Sample Renyi function with the fitted Model 6 for a large window of sky data.
plot(
  Tq[, 1],
  Tq[, 2],
  ylab = "T(q)",
  xlab = "q",
  pch = 20,
  col = "blue")
lines(Tq[, 1], Q16, col = "red", lwd = 2)

########### Figure 8e-Difference between the Sample Renyi function and the fitted Model 6
#Figure 8e gives the plot of the difference between the Sample Renyi function and the fitted Model 6 for a large window of sky data.
plot(
  Tq[, 1],
  Q16 - Tq[, 2],
  ylab = "Difference",
  xlab = "q",
  pch = 20,
  col = "blue")

residuals <- Q16 - Tq[, 2]
sqrt(mean(residuals ^ 2))

