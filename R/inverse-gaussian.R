## --##--##--##--##--##--##--##--##--##--##--##--##--##--##--##
##   Functions related to the inverse gaussian distribution  ##
## --##--##--##--##--##--##--##--##--##--##--##--##--##--##--##

#' @param m vector of means
#' @param s vector of dispersion parameters
#' @rdname rtrunc
#' @export
rtruncinvgauss <- rtrunc.invgauss <- function(n, m, s, a = 0, b = Inf) {
  class(n) <- "trunc_invgauss"
  sampleFromTruncated(mget(ls()))
}

#' @export
dtrunc.trunc_invgauss <- function(y, m, s, eta, a = 0, b = Inf, ...) {
  if (missing(eta)) {
    eta <- parameters2natural.parms_invgauss(c("m" = m, "s" = s))
  }
  parm <- natural2parameters.parms_invgauss(eta)
  dens <- rescaledDensities(y, a, b, dinvgauss, pinvgauss, parm["m"], parm["s"])
  return(dens)
}

#' @rdname dtrunc
#' @export
dtruncinvgauss <- dtrunc.trunc_invgauss

#' @export
empiricalParameters.trunc_invgauss <- function(y, ...) {
  # Returns empirical parameter estimates mean and shape
  mean <- mean(y)
  sd <- sd(y)
  lambda <- mean ^ 3 / sd ^ 2
  parms <- c(m = mean, s = 1 / lambda)
  class(parms) <- "parms_invgauss"
  return(parms)
}

sufficientT.trunc_invgauss <- function(y) {
  return(suff.T = cbind(y, 1 / y))
}

#' @export
parameters2natural.parms_invgauss <- function(parms, ...) {
  # parms: The parameters mean and shape in a normal distribution
  # returns the natural parameters
  mu <- parms[["m"]]
  lambda <- 1 / parms[["s"]]
  eta <- c(eta1 = -lambda / (2 * mu ^ 2), eta2 = -lambda / 2)
  class(eta) <- class(parms)
  return(eta)
}

#' @export
natural2parameters.parms_invgauss <- function(eta, ...) {
  # eta: The natural parameters in an inverse gaussian distribution
  # returns (mean,shape)
  if (length(eta) != 2) stop("Eta must be a vector of two elements")
  mu <- sqrt(eta[[2]] / eta[[1]])
  lambda <- -2 * eta[[2]]
  parms <- c(m = mu, s = 1 / lambda)
  class(parms) <- class(eta)
  return(parms)
}

getYseq.trunc_invgauss <- function(y, y.min, y.max, n = 100) {
  m <- mean(y, na.rm = TRUE)
  sd <- sd(y, na.rm = TRUE)
  lo <- max(0, y.min, m - 3.5 * sd)
  hi <- min(y.max, m + 3.5 * sd)
  out <- seq(lo, hi, length = n)
  out <- out[out > 0] # y must be positive
  class(out) <- class(y)
  return(out)
}

getGradETinv.parms_invgauss <- function(eta, ...) {
  # eta: Natural parameter
  # return the inverse of E.T differentiated with respect to eta' : p x p matrix
  mx_11 <- -sqrt(eta[2] / eta[1] ^ 3)
  mx_12 <- 1 / (sqrt(eta[1] * eta[2]))
  mx_21 <- mx_12
  mx_22 <- (1 - eta[1] * sqrt(eta[2] / eta[1])) / (eta[2] ^ 2)
  A_inv <- 0.5 * matrix(c(mx_11, mx_12, mx_21, mx_22), ncol = 2)
  return(A = solve(A_inv))
}
