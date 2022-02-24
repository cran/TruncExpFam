## --##--##--##--##--##--##--##--##--##--##--##--##--##
##   Functions related to the Beta distribution      ##
##         Variant 1                                 ##
## --##--##--##--##--##--##--##--##--##--##--##--##--##

#' @importFrom stats rbeta
#' @param shape1 positive shape parameter alpha
#' @param shape2 positive shape parameter beta
#' @rdname rtrunc
#' @export
rtruncbeta <- rtrunc.beta <- function(n, shape1, shape2, a = 0, b = 1) {
  class(n) <- "trunc_beta"
  sampleFromTruncated(mget(ls()))
}

#' @export
dtrunc.trunc_beta <- function(y, eta, a = 0, b = 1) {
  parm <- natural2parameters.trunc_beta(eta)
  dens <- ifelse((y < a) | (y > b), 0, dbeta(y, shape1 = parm[1], shape2 = parm[2]))
  if (!missing(a)) {
    F.a <- pbeta(a, shape1 = parm[1], shape2 = parm[2])
  } else {
    F.a <- 0
  }
  if (!missing(b)) {
    F.b <- pbeta(b, shape1 = parm[1], shape2 = parm[2])
  } else {
    F.b <- 1
  }
  const <- 1 / (F.b - F.a)
  return(dens * const)
}

#' @importFrom stats dbeta pbeta
#' @rdname dtrunc
#' @export
dtruncbeta <- dtrunc.trunc_beta

#' @export
init.parms.trunc_beta <- function(y) {
  # Returns  parameter estimates mean and sd
  amean <- mean(y)
  avar <- var(y)
  alpha <- amean^2 * (1 - amean) / avar - amean
  beta <- alpha * (1 / amean - 1)
  parms <- c(shape1 = alpha, shape2 = beta)
  class(parms) <- "trunc_beta"
  return(parms)
}

sufficientT.trunc_beta <- function(y) {
  return(suff.T = cbind(log(y), log(1 - y)))
}

averageT.trunc_beta <- function(y) {
  return(apply(cbind(log(y), log(1 - y)), 2, mean))
}

#' @export
natural2parameters.trunc_beta <- function(eta) {
  # eta: The natural parameters in a beta distribution
  # returns (alpha,beta)
  parms <- c(shape1 = eta[1], shape2 = eta[2])
  class(parms) <- class(eta)
  return(parms)
}

#' @export
parameters2natural.trunc_beta <- function(parms) {
  # parms: The parameters shape and rate in a beta distribution
  # returns the natural parameters
  eta <- c(shape1 = parms[1], shape2 = parms[2])
  class(eta) <- class(parms)
  return(eta)
}

getYseq.trunc_beta <- function(y, y.min = 0, y.max = 1, n = 100) {
  # needs chekking
  mean <- mean(y, na.rm = TRUE)
  sd <- var(y, na.rm = TRUE)^0.5
  lo <- max(y.min, mean - 5 * sd)
  hi <- min(y.max, mean + 5 * sd)
  out <- seq(lo, hi, length = n)
  class(out) <- class(y)
  return(out)
}

getGradETinv.trunc_beta <- function(eta) {
  # eta: Natural parameter
  # return the inverse of E.T differentiated with respect to eta' : p x p matrix
  term.1 <- sum(1 / (((1:10000) + eta[1]))^2)
  term.2 <- sum(1 / (((1:10000) + eta[2]))^2)
  term.12 <- sum(1 / (((1:10000) + eta[1] + eta[2]))^2)
  return(A = solve(matrix(c(term.1 - term.12, -term.12, -term.12, term.2 - term.12, ncol = 2))))
}
