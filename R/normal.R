## --##--##--##--##--##--##--##--##--##--##--##--##--##
##   Functions related to the normal distribution   ##
## --##--##--##--##--##--##--##--##--##--##--##--##--##

#' @param mean mean of parent distribution
#' @param sd standard deviation is parent distribution
#' @rdname rtrunc
#' @export
rtruncnorm <- rtrunc.normal <- function(n, mean, sd, a = -Inf, b = Inf) {
  class(n) <- "trunc_normal"
  sampleFromTruncated(mget(ls()))
}

#' @export
dtrunc.trunc_normal <- function(y, eta, a = -Inf, b = Inf) {
  parm <- natural2parameters.trunc_normal(eta)
  dens <- ifelse((y < a) | (y > b), 0, dnorm(y, mean = parm[1], sd = parm[2]))
  if (!missing(a)) {
    F.a <- pnorm(a, parm[1], parm[2])
  } else {
    F.a <- 0
  }
  if (!missing(b)) {
    F.b <- pnorm(b, parm[1], parm[2])
  } else {
    F.b <- 1
  }
  const <- 1 / (F.b - F.a)
  return(dens * const)
}

#' @rdname dtrunc
#' @export
dtruncnorm <- dtrunc.trunc_normal

#' @export
init.parms.trunc_normal <- function(y) {
  # Returns empirical parameter estimates mean and sd
  parms <- c(mean = mean(y), sd = sqrt(var(y)))
  class(parms) <- "trunc_normal"
  return(parms)
}

sufficientT.trunc_normal <- function(y) {
  return(suff.T = cbind(y, y^2))
}

averageT.trunc_normal <- function(y) {
  return(apply(sufficientT.trunc_normal(y), 2, mean))
}

#' @export
natural2parameters.trunc_normal <- function(eta) {
  # eta: The natural parameters in a normal distribution
  # returns (mean,sigma)
  parms <- c(mean = -0.5 * eta[1] / eta[2], sd = sqrt(-0.5 / eta[2]))
  class(parms) <- class(eta)
  return(parms)
}

#' @export
parameters2natural.trunc_normal <- function(parms) {
  # parms: The parameters mean and sd in a normal distribution
  # returns the natural parameters
  eta <- c(eta.1 = parms[1], eta.2 = -0.5) / parms[2]^2
  class(eta) <- class(parms)
  return(eta)
}

getYseq.trunc_normal <- function(y, y.min, y.max, n = 100) {
  mean <- mean(y, na.rm = TRUE)
  sd <- var(y, na.rm = TRUE)^0.5
  lo <- max(y.min, mean - 3.5 * sd)
  hi <- min(y.max, mean + 3.5 * sd)
  out <- seq(lo, hi, length = n)
  class(out) <- class(y)
  return(out)
}

getGradETinv.trunc_normal <- function(eta, ...) {
  # eta: Natural parameter
  # return the inverse of E.T differentiated with respect to eta' : p x p matrix
  return(A = solve(0.5 * matrix(c(-1 / eta[2], eta[1] / eta[2]^2, eta[1] / eta[2]^2, 1 / eta[2]^2 - eta[1]^2 / eta[2]^3), ncol = 2)))
}
