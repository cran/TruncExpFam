## --##--##--##--##--##--##--##--##--##--##--##--##--##
##   Functions related to the Poisson distribution   ##
## --##--##--##--##--##--##--##--##--##--##--##--##--##

#' @param lambda mean and var of "parent" distribution
#' @rdname rtrunc
#' @export
rtruncpois <- rtrunc.poisson <- function(n, lambda, a = 0, b = Inf) {
  class(n) <- "trunc_poisson"
  sampleFromTruncated(mget(ls()))
}

#' @export
dtrunc.trunc_poisson <- function(y, eta, a = 0, b = Inf) {
  parm <- exp(eta)
  dens <- ifelse((y < a) | (y > b), 0, dpois(y, parm))
  if (!missing(a)) {
    F.a <- ppois(a - 1, parm)
  } else {
    F.a <- 0
  }
  if (!missing(b)) {
    F.b <- ppois(b, parm)
  } else {
    F.b <- 1
  }
  return(dens / (F.b - F.a))
}

#' @rdname dtrunc
#' @export
dtruncpois <- dtrunc.trunc_poisson

#' @export
init.parms.trunc_poisson <- function(y) {
  # Returns empirical parameter estimate for lambda
  parms <- mean(y)
  class(parms) <- "trunc_poisson"
  return(parms)
}

sufficientT.trunc_poisson <- function(y) {
  return(suff.T = y)
}

averageT.trunc_poisson <- function(y) {
  return(mean(y))
}


#' @export
natural2parameters.trunc_poisson <- function(eta) {
  # eta: The natural parameters in a Poisson distribution
  # returns (mean,sigma)
  lambda <- c(lambda = exp(eta))
  class(lambda) <- class(eta)
  return(lambda)
}

#' @export
parameters2natural.trunc_poisson <- function(parms) {
  # parms: The parameter lambda in a Poisson distribution
  # returns the natural parameters
  eta <- log(parms)
  class(eta) <- class(parms)
  return(eta)
}

getGradETinv.trunc_poisson <- function(eta, ...) {
  # eta: Natural parameter
  # return the inverse of E.T differentiated with respect to eta
  return(A = exp(-eta))
}

getYseq.trunc_poisson <- function(y, y.min = 0, y.max, n = 100) {
  mean <- mean(y, na.rm = TRUE)
  var.y <- var(y, na.rm = TRUE)
  lo <- max(round(y.min), 0)
  hi <- min(y.max, round(mean + 10 * sqrt(var.y)))
  out <- seq(lo, hi)
  class(out) <- class(y)
  return(out)
}
