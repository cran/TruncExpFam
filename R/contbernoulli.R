## --##--##--##--##--##--##--##--##--##--##--##--##--##--##--##--##
##   Functions related to the continuous Bernoulli distribution  ##
## --##--##--##--##--##--##--##--##--##--##--##--##--##--##--##--##

# Sampling function for a continuous bernoulli distribution
# This distribution is not implemented in Base R
# Used in the sampling of the truncated continuous bernoulli
rcontbern <- function(n, lambda) {
  if ((lambda < 0) || (lambda > 1)) {
    stop("lambda must be in (0, 1)")
  }
  u <- runif(n)
  if (lambda == 0.5) {
    return(u)
  }

  # The inverse of the CDF for a cont. bernoulli distribution
  x <- log(1 + (2 * lambda - 1) * u / (1 - lambda)) /
    log(lambda / (1 - lambda))
  return(x)
}

#' @param lambda mean of "parent" distribution
#' @rdname rtrunc
#' @export
rtrunccontbern <- rtrunc.contbern <- function(n, lambda, a = 0, b = 1) {
  class(n) <- "trunc_contbern"
  sampleFromTruncated(mget(ls()))
}

# The two functions 'dcontbern' and 'pcontbern' below act in support of the
# truncated continuous bernoulli distribution, as base R does not include
# this distribution
# dcontbern is the untruncated function (which is not present in base R)
dcontbern <- function(x, lambda) {
  norm.const <- ifelse(
    test = lambda == 0.5,
    yes  = 2,
    no   = 2 * (atanh(1 - 2 * lambda)) / (1 - 2 * lambda)
  )
  d <- norm.const * (lambda^x) * (1 - lambda) ^ (1 - x)
  class(d) <- class(x)
  return(d)
}

qcontbern <- function(p, lambda) {
  if (lambda == .5) {
    return(p)
  } else {
    term1 <- log(2 * lambda * p - p + 1 - lambda)
    term2 <- log(1 - lambda)
    term3 <- log(lambda)
    return((term1 - term2) / (term3 - term2))
  }
}

# untruncated version (not implemented in base R)
pcontbern <- function(x, lambda) {
  p <- ((lambda^x) * (1 - lambda) ^ (1 - x) + lambda - 1) / (2 * lambda - 1)
  return(p)
}

#' @export
dtrunc.trunc_contbern <- function(
  y, lambda, eta, a = 0, b = 1, ...
) {
  if (missing(eta)) {
    eta <- parameters2natural.parms_contbern(c("lambda" = lambda))
  }
  lambda <- natural2parameters.parms_contbern(eta)
  dens <- rescaledDensities(y, a, b, dcontbern, pcontbern, lambda)
  return(dens)
}

#' @rdname dtrunc
#' @export
dtrunccontbern <- dtrunc.trunc_contbern

#' @export
#' @param eta vector of natural parameters
#' @rdname dtrunc
dtrunccontbern <- dtrunc.trunc_contbern

#' @export
empiricalParameters.trunc_contbern <- function(y, ...) {
  # Returns empirical parameter estimate for the lambda parameter
  # Note: lambda cannot be expressed in closed form as a function of the mean
  parms <- c("lambda" = mean(y))
  class(parms) <- "parms_contbern"
  return(parms)
}

sufficientT.trunc_contbern <- function(y) {
  return(suff.T = y)
}

#' @export
natural2parameters.parms_contbern <- function(eta, ...) {
  # eta: The natural parameters in a continuous bernoulli distribution
  # returns rate
  if (length(eta) != 1) stop("Eta must be one single number")
  rate <- c(lambda = 1 / (1 + exp(-eta[[1]])))
  class(rate) <- class(eta)
  return(rate)
}

#' @export
parameters2natural.parms_contbern <- function(parms, ...) {
  # parms: The parameter lambda in a continuous bernoulli distribution
  # returns the natural parameters
  eta <- prepEta(log(parms / (1 - parms)), class(parms))
  return(eta)
}

getYseq.trunc_contbern <- function(y, y.min = 0, y.max, n = 100) {
  mean <- mean(y, na.rm = TRUE)
  var.y <- var(y, na.rm = TRUE)
  lo <- max(round(y.min), 0)
  hi <- min(y.max, round(mean + 10 * sqrt(var.y)), 1)
  out <- seq(lo, hi, length = n)
  class(out) <- class(y)
  return(out)
}

getGradETinv.parms_contbern <- function(eta, ...) {
  # eta: Natural parameter
  # return the inverse of E.T differentiated with respect to eta
  exp.eta <- exp(eta)
  return(A = ((exp.eta - 1) * eta)^2 / (exp.eta * (exp.eta - eta^2 + eta - 1)))
}
