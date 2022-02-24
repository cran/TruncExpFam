## --##--##--##--##--##--##--##--##--##--##--##--##--##--##--##--##
##   Functions related to the continuous Bernoulli distribution  ##
## --##--##--##--##--##--##--##--##--##--##--##--##--##--##--##--##

#' @importFrom stats runif
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
  if ((x < 0) | (x > 1)) {
    return(0)
  }
  if (lambda == 0.5) {
    norm.const <- 2
  } else {
    norm.const <- 2 * (atanh(1 - 2 * lambda)) / (1 - 2 * lambda)
  }
  d <- norm.const * (lambda^x) * (1 - lambda) ^ (1 - x)
  return(d)
}

# untruncated version (not implemented in base R)
pcontbern <- function(x, lambda) {
  if (x < 0) {
    p <- 0
  } else if (x > 1) {
    p <- 1
  } else if (lambda == 0.5) {
    p <- x
  } else {
    p <- ((lambda^x) * (1 - lambda) ^ (1 - x) + lambda - 1) / (2 * lambda - 1)
  }
  return(p)
}

#' @rdname dtrunc
#' @export
dtrunccontbern <- dtrunc.trunc_contbern <- function(y, eta, a = 0, b = 1) {
  lambda <- natural2parameters.trunc_contbern(eta)
  dens <- ifelse((y <= a) | (y > b), 0, dcontbern(y, lambda = lambda))
  if (!missing(a)) {
    F.a <- pcontbern(a, lambda)
  } else {
    F.a <- 0
  }
  if (!missing(b)) {
    F.b <- pcontbern(b, lambda)
  } else {
    F.b <- 1
  }
  return(dens / (F.b - F.a))
}

#' @export
init.parms.trunc_contbern <- function(y) {
  # Returns empirical parameter estimate for the lambda parameter
  # Note: lambda cannot be expressed in closed form as a function of the mean
  parms <- mean(y)
  class(parms) <- "trunc_contbern"
  return(parms)
}

sufficientT.trunc_contbern <- function(y) {
  return(suff.T = y)
}

averageT.trunc_contbern <- function(y) {
  return(mean(y))
}

#' @export
natural2parameters.trunc_contbern <- function(eta) {
  # eta: The natural parameters in a continuous bernoulli distribution distribution
  # returns rate
  rate <- c(lamda = 1 / (1 + exp(-eta)))
  class(rate) <- class(eta)
  return(rate)
}

#' @export
parameters2natural.trunc_contbern <- function(parms) {
  # parms: The parameter lambda in a continuous bernoulli distribution
  # returns the natural parameters
  eta <- log(parms / (1 - parms))
  class(eta) <- class(parms)
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

getGradETinv.trunc_contbern <- function(eta) {
  # eta: Natural parameter
  # return the inverse of E.T differentiated with respect to eta
  exp.eta <- exp(eta)
  return(A = ((exp.eta - 1) * eta)^2 / (exp.eta * (exp.eta - eta^2 + eta - 1)))
}
