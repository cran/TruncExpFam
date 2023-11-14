## --##--##--##--##--##--##--##--##--##--##--##--##--##
##   Functions related to the Beta distribution      ##
##         Variant 1                                 ##
## --##--##--##--##--##--##--##--##--##--##--##--##--##

#' @param shape1 positive shape parameter alpha
#' @param shape2 positive shape parameter beta
#' @rdname rtrunc
#' @export
rtrunc.beta <- function(n, shape1, shape2, a = 0, b = 1) {
  class(n) <- "trunc_beta"
  sampleFromTruncated(mget(ls()))
}

#' @rdname rtrunc
#' @export
rtruncbeta <- rtrunc.beta

#' @export
dtrunc.trunc_beta <- function(y, shape1, shape2, eta, a = 0, b = 1, ...) {
  if (missing(eta)) {
    eta <- parameters2natural.parms_beta(c(shape1, shape2))
  }
  parm <- natural2parameters.parms_beta(eta)
  dens <- rescaledDensities(y, a, b, dbeta, pbeta, parm[1], parm[2])
  return(dens)
}

#' @rdname dtrunc
#' @inheritParams rtrunc
#' @export
dtruncbeta <- dtrunc.trunc_beta

#' @export
empiricalParameters.trunc_beta <- function(y, ...) {
  # Returns  parameter estimates mean and sd
  amean <- mean(y)
  avar <- var(y)
  alpha <- amean^2 * (1 - amean) / avar - amean
  beta <- alpha * (1 / amean - 1)
  parms <- c(shape1 = alpha, shape2 = beta)
  class(parms) <- "parms_beta"
  return(parms)
}

sufficientT.trunc_beta <- function(y) {
  # Calculates the sufficient statistic T(y)
  return(suff.T = cbind(log(y), log(1 - y)))
}

#' @export
natural2parameters.parms_beta <- function(eta, ...) {
  # eta: The natural parameters in a beta distribution
  # returns (alpha,beta)
  if (length(eta) != 2) stop("Eta must be a vector of two elements")
  parms <- c(shape1 = eta[[1]], shape2 = eta[[2]])
  class(parms) <- class(eta)
  return(parms)
}

#' @export
parameters2natural.parms_beta <- function(parms, ...) {
  # parms: The parameters shape and rate in a beta distribution
  # returns the natural parameters
  eta <- prepEta(c(parms[1], parms[2]), class(parms))
  return(eta)
}

getYseq.trunc_beta <- function(y, y.min = 0, y.max = 1, n = 100) {
  # needs chekking
  mean <- mean(y, na.rm = TRUE)
  sd <- var(y, na.rm = TRUE)^0.5
  lo <- max(y.min, mean - 5 * sd)
  hi <- min(y.max, mean + 5 * sd)
  out <- seq(lo, hi, length = n)
  out <- out[out > 0 & out < 1] # prevents NaN as sufficient statistics
  class(out) <- class(y)
  return(out)
}

getGradETinv.parms_beta <- function(eta, ...) {
  # eta: Natural parameter
  # return the inverse of E.T differentiated with respect to eta' : p x p matrix
  # Uses approximation for the digamma function: digamma(x) ~ ln(x) - 1 / 2 / x
  # Source: https://en.wikipedia.org/wiki/Digamma_function
  # Derivatives with respect to etas calculated on wolfram alpha
  x <- eta[1]
  y <- eta[2]
  term.1 <- (y * (2 * x ^ 2 + y + 2 * x * (1 + y))) / (2 * x ^ 2 * (x + y) ^ 2)
  term.12 <- -(1 + 2 * (x + y)) / (2 * (x + y) ^ 2)
  term.2 <- (x * (x + 2 * x * y + 2 * y * (1 + y))) / (2 * y ^ 2 * (x + y) ^ 2)
  A_inv <- matrix(c(term.1, term.12, term.12, term.2), ncol = 2)
  return(A = solve(A_inv))
}
