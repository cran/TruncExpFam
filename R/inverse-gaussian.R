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
dtrunc.trunc_invgauss <- function(y, eta, a = 0, b = Inf) {
  parm <- natural2parameters.trunc_invgauss(eta)
  dens <- ifelse((y < a) | (y > b), 0, dinvgauss(y, m = parm[1], s = parm[2]))
  if (!missing(a)) {
    # ERROR		F.a <- pnorm(a, parm[1], parm[2])
    F.a <- pinvgauss(a, parm[1], parm[2])
  } else {
    F.a <- 0
  }
  if (!missing(b)) {
    # ERROR		F.b <- pnorm(b, parm[1], parm[2])
    F.b <- pinvgauss(b, parm[1], parm[2])
  } else {
    F.b <- 1
  }
  const <- 1 / (F.b - F.a)
  return(dens * const)
}

#' @importFrom rmutil rinvgauss
#' @rdname dtrunc
#' @export
#' @importFrom rmutil dinvgauss pinvgauss
dtruncinvgauss <- dtrunc.trunc_invgauss

#' @export
init.parms.trunc_invgauss <- function(y) {
  # Returns empirical parameter estimates mean and shape
  mean <- mean(y)
  shp <- 1 / (mean(1 / y) - 1 / mean)
  parms <- c(mean = mean, shape = shp)
  class(parms) <- "trunc_invgauss"
  return(parms)
}

sufficientT.trunc_invgauss <- function(y) {
  return(suff.T = cbind(y, 1 / y))
}

averageT.trunc_invgauss <- function(y) {
  return(apply(sufficientT.trunc_invgauss(y), 2, mean))
}

#' @export
natural2parameters.trunc_invgauss <- function(eta) {
  # eta: The natural parameters in an inverse gaussian distribution
  # returns (mean,shape)
  parms <- c(mean = sqrt(eta[2] / eta[1]), shape = -2 * eta[2])
  class(parms) <- class(eta)
  return(parms)
}

#' @export
parameters2natural.trunc_invgauss <- function(parms) {
  # parms: The parameters mean and shape in a normal distribution
  # returns the natural parameters
  eta <- c(eta.1 = -parms[2] / (2 * parms[1]^2), eta.2 = -0.5 * parms[2])
  class(eta) <- class(parms)
  return(eta)
}

getYseq.trunc_invgauss <- function(y, y.min, y.max, n = 100) {
  mean <- mean(y, na.rm = TRUE)
  shape <- var(y, na.rm = TRUE)^0.5
  lo <- max(max(0, y.min), mean - 3.5 * shape)
  hi <- min(y.max, mean + 3.5 * shape)
  out <- seq(lo, hi, length = n)
  class(out) <- class(y)
  return(out)
}

getGradETinv.trunc_invgauss <- function(eta, ...) {
  # eta: Natural parameter
  # return the inverse of E.T differentiated with respect to eta' : p x p matrix
  sqrt.eta1 <- sqrt(eta[1])
  sqrt.eta2 <- sqrt(eta[2])
  return(A = solve(0.5 * matrix(c(
    -sqrt.eta2 / sqrt.eta1^3, 1 / (sqrt.eta1 * sqrt.eta2), 1 / (sqrt.eta1 * sqrt.eta2),
    -sqrt.eta1 / sqrt.eta2^3 + 1 / sqrt.eta2^2
  ), ncol = 2)))
}
