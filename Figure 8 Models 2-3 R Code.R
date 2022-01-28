setwd(dirname(rstudioapi::getActiveDocumentContext()$path))

library(stats)
library(rcosmo)
library(sp)
library(minpack.lm)

############## Renyi function
#The function fRen computes the Renyi function for a chosen window of sky data(e.g.:- whole sky, large, medium windows,etc.).
fRen1 <-
  function(cmbdf,
           q.min = 1.01,
           q.max = 10,
           N = 20,
           k.box = log2(nside(cmbdf)) - 3,
           intensities = "I"){
    if (!is.CMBDataFrame(cmbdf)) {
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
    if (res.max - k.box > 0) {
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

############### Model 2: Large Window near the pole of the sphere  #############

cmbdf <- CMBDataFrame("CMB_map_smica1024.fits")
minint <- min(cmbdf[, "I", drop = T])
q.min <- 1.01
q.max <- 2.0
N <- 20
#Here, a large window of CMB sky data is chosen such that the area=1.231.
win <-
  CMBWindow(theta = c(3 * pi / 6, 3 * pi / 6, pi / 4, pi / 4),
            phi = c(0, pi / 2, pi / 2, 0))
cmbdf1 <- window(cmbdf, new.window = win)

Tq <- fRen1(cmbdf, q.min, q.max, N)
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

QM2 <-
  nlsLM(
    y ~  (log(1 - (x * B)) - x * log(1 - B))/A,
    start = list(A = 2, B = 0.1),
    data = Ren1,
    trace = TRUE,
    control = nls.control(
      maxiter = 1000,
      tol = 1e-2,
      printEval = FALSE,
      warnOnly = TRUE))
#Coef(QM2) gives the estimated parameters A and B
coef(QM2)

#Verifying that the estimated parameters satisfy the assumptions of Theorem 6.2
#Evaluating b(L.H.S)
exp(coef(QM2)[1])
#Evaluating R.H.S
1+coef(QM2)[2]^2/(1-2*coef(QM2)[2])

Q21 <- (fitted(QM2) + x - 1)

#This figure gives the plot of the sample Renyi function with the fitted log-gamma model for a large window of sky data
plot(
  Tq[, 1],
  Tq[, 2],
  ylab = "T(q)",
  xlab = "q",
  pch = 20,
  col = "blue")
lines(Tq[, 1], Q21, col = "red", lwd = 2)

########### Figure 8a-Difference between the Sample Renyi function and the fitted log-gamma model
#Figure 8a gives the plot of the difference between the Sample Renyi function and the fitted log-gamma model for a large window of sky data.
plot(
  Tq[, 1],
  Q21 - Tq[, 2],
  ylab = "Difference",
  xlab = "q",
  pch = 20,
  col = "blue")

residuals <- Q21 - Tq[, 2]
sqrt(mean(residuals ^ 2))

############# Model 3: Large Window near the pole of the sphere ##########

QM3 <-
  nlsLM(
    y ~ 1/A * (
      x * log(2*(C^B)*besselK(2 * C, B, expon.scaled = FALSE) / gamma(B)) - B * log(x) / 2 - log(besselK(2 * C * sqrt(x), B, expon.scaled = FALSE)
      ) - log(2 *(C^ B) / gamma(B))),
    start = list(A = 0.1,  B = 2, C =1.4),
    data = Ren1,
    trace = TRUE,
    lower    = c(0.01, 0.1, 0.1),
    upper =    c(10, 10, 10),
    control = nls.control(maxiter =1000, tol = 1e-2,
      printEval = FALSE))
#Coef(QM3) gives the estimated parameters A, B and C.
coef(QM3)

#Verifying that the estimated parameters satisfy the assumptions of Theorem 6.3
#Evaluating b(L.H.S)
exp((coef(QM3)[1])/2)
#Evaluating R.H.S
sqrt((gamma(coef(QM3)[2])*2^((coef(QM3)[2]/2)-1)*besselK(2*sqrt(2)*(coef(QM3)[3]), coef(QM3)[2], expon.scaled = FALSE))/(((coef(QM3)[3])^coef(QM3)[2])*((besselK((2*(coef(QM3)[3])), coef(QM3)[2], expon.scaled = FALSE))^2)))

Q31 <- (fitted(QM3) + x - 1)

#This figure gives the plot of the sample Renyi function with the fitted log-negative-inverse-gamma model for a large window of sky data.
plot(
  Tq[, 1],
  Tq[, 2],
  ylab = "T(q)",
  xlab = "q",
  pch = 20,
  col = "blue")
lines(Tq[, 1], Q31, col = "red", lwd = 2)

########### Figure 8b-Difference between the Sample Renyi function and the fitted log-negative-inverse-gamma model
#Figure 8b gives the plot of the difference between the Sample Renyi function and the fitted log-negative-inverse-gamma model for a large window of sky data.
plot(
  Tq[, 1],
  Q31 - Tq[, 2],
  ylab = "Difference",
  xlab = "q",
  pch = 20,
  col = "blue")

residuals <- Q31 - Tq[, 2]
sqrt(mean(residuals ^ 2))

