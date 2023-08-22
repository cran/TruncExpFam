## --##--##--##--##--##--##--##--##--##--##--##--##--##
##   Functions related to the Binomial distribution  ##
## --##--##--##--##--##--##--##--##--##--##--##--##--##

#' @param size number of trials
#' @param prob probability of success on each trial
#' @rdname rtrunc
#' @export
rtruncbinom <- rtrunc.binomial <- function(n, size, prob, a = 0, b = size) {
  class(n) <- "trunc_binomial"
  sampleFromTruncated(mget(ls()))
}

#' @export
dtrunc.trunc_binomial <- function(y, eta, a = 0, b = attr(y, "parameters")$size, ...) {
  nsize <- attr(y, "parameters")$size
  my.dbinom <- function(nsize) dbinom(y, size = nsize, prob = proba)
  my.pbinom <- function(z, nsize) pbinom(z, size = nsize, prob = proba)
  proba <- 1 / (1 + exp(-eta))
  dens <- ifelse((y < a) | (y > b), 0, my.dbinom(nsize))
  if (!missing(a)) {
    F.a <- my.pbinom(a - 1, nsize)
  } else {
    F.a <- 0
  }
  if (!missing(b)) {
    F.b <- my.pbinom(b, nsize)
  } else {
    F.b <- 1
  }
  return(dens / (F.b - F.a))
}

#' @rdname dtrunc
#' @param ... size
#' @export
dtruncbinom <- dtrunc.trunc_binomial

#' @export
init.parms.trunc_binomial <- function(y) {
  # Returns empirical parameter estimate for lambda
  nsize <- attr(y, "parameters")$size
  parms <- mean(y / nsize)
  attr(parms, "nsize") <- nsize
  class(parms) <- "trunc_binomial"
  return(parms)
}

sufficientT.trunc_binomial <- function(y) {
  return(suff.T = y)
}

averageT.trunc_binomial <- function(y) {
  return(mean(y))
}

density.trunc_binomial <- function(y, eta, ...) {
  parms <- 1 / (1 + exp(-eta))
  dbinom(y, ..., parms)
}

#' @export
natural2parameters.trunc_binomial <- function(eta) {
  # eta: The natural parameters in a binomial distribution
  # returns (p)
  p <- c(prob = 1 / (1 + exp(-eta)))
  class(p) <- class(eta)
  return(p)
}

#' @export
parameters2natural.trunc_binomial <- function(parms) {
  # parms: The probability parameter p in a binomial distribution
  # returns the natural parameters
  eta <- log(parms / (1 - parms))
  class(eta) <- class(parms)
  return(eta)
}

getGradETinv.trunc_binomial <- function(eta, ...) {
  # eta: Natural parameter
  # return the inverse of E.T differentiated with respect to eta
  nsize <- attr(eta, "nsize")
  exp.eta <- exp(eta)
  return(A = ((1 + exp.eta)^2 / exp.eta) / nsize)
}

getYseq.trunc_binomial <- function(y, y.min = 0, y.max, n = 100) {
  nsize <- attr(y, "parameters")$size
  y.lo <- round(y.min)
  y.hi <- round(y.max)
  lo <- max(y.lo, 0)
  hi <- min(y.hi, nsize)
  out <- seq(lo, hi)
  attributes(out) <- attributes(y)
  return(out)
}
