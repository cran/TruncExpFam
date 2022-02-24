## --##--##--##--##--##--##--##--##--##--##--##--##--##--##
##   Functions related to the Chi Square distribution    ##
## --##--##--##--##--##--##--##--##--##--##--##--##--##--##

#' @importFrom stats rexp
#' @param rate vector of rates
#' @rdname rtrunc
#' @export
rtruncexp <- rtrunc.exp <- function(n, rate = 1, a = 0, b = Inf) {
  class(n) <- "trunc_exp"
  sampleFromTruncated(mget(ls()))
}

#' @export
dtrunc.trunc_exp <- function(y, eta, a = 0, b = Inf) {
  rate <- natural2parameters.trunc_exp(eta)
  dens <- ifelse((y <= a) | (y > b), 0, dexp(y, rate = rate))
  if (!missing(a)) {
    F.a <- pexp(a, rate)
  } else {
    F.a <- 0
  }
  if (!missing(b)) {
    F.b <- pexp(b, rate)
  } else {
    F.b <- 1
  }
  return(dens / (F.b - F.a))
}

#' @importFrom stats dexp pexp
#' @rdname dtrunc
#' @export
dtruncexp <- dtrunc.trunc_exp

#' @export
init.parms.trunc_exp <- function(y) {
  # Returns empirical parameter estimate for the rate parameter
  parms <- mean(y)
  class(parms) <- "trunc_exp"
  return(parms)
}

sufficientT.trunc_exp <- function(y) {
  return(suff.T = y)
}

averageT.trunc_exp <- function(y) {
  return(mean(y))
}

#' @export
natural2parameters.trunc_exp <- function(eta) {
  # eta: The natural parameters in an exponential distribution distribution
  # returns rate
  lambda <- -eta
  class(lambda) <- class(eta)
  return(lambda)
}

#' @export
parameters2natural.trunc_exp <- function(parms) {
  # parms: The parameter lambda in an exponential distribution
  # returns the natural parameters
  eta <- -parms
  class(eta) <- class(parms)
  return(eta)
}

getGradETinv.trunc_exp <- function(eta) {
  # eta: Natural parameter
  # return the inverse of E.T differentiated with respect to eta
  return(A = eta^2)
}

getYseq.trunc_exp <- function(y, y.min = 0, y.max, n = 100) {
  mean <- mean(y, na.rm = TRUE)
  var.y <- var(y, na.rm = TRUE)
  lo <- max(round(y.min), 0)
  hi <- min(y.max, round(mean + 10 * sqrt(var.y)))
  out <- seq(lo, hi, length = n)
  class(out) <- class(y)
  return(out)
}
