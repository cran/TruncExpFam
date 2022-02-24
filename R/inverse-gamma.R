## --##--##--##--##--##--##--##--##--##--##--##--##--##--##--##
##   Functions related to the inverse gamma distribution     ##
## --##--##--##--##--##--##--##--##--##--##--##--##--##--##--##

#' @importFrom invgamma rinvgamma
#' @param shape inverse gamma shape parameter
#' @param rate inverse gamma rate parameter
#' @param scale inverse gamma scale parameter
#' @rdname rtrunc
#' @export
rtruncinvgamma <- rtrunc.invgamma <- function(n, shape, rate = 1, scale = 1 / rate,
                                              a = 0, b = Inf) {
  class(n) <- "trunc_invgamma"
  sampleFromTruncated(mget(ls()))
}

#' @importFrom invgamma dinvgamma pinvgamma
#' @export
dtrunc.trunc_invgamma <- function(y, eta, a = 0, b = Inf) {
  parm <- natural2parameters.trunc_invgamma(eta)
  dens <- ifelse((y < a) | (y > b), 0, dinvgamma(y, shape = parm[1], rate = parm[2]))
  if (!missing(a)) {
    F.a <- pinvgamma(a, shape = parm[1], rate = parm[2])
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

#' @rdname dtrunc
#' @export
dtruncinvgamma <- dtrunc.trunc_invgamma

#' @export
init.parms.trunc_invgamma <- function(y) {
  # Returns  parameter estimates mean and sd
  amean <- mean(y)
  avar <- var(y)
  alpha <- amean^2 / avar + 2
  beta <- (alpha - 1) * amean
  parms <- c(shape = alpha, rate = beta)
  class(parms) <- "trunc_invgamma"
  return(parms)
}

sufficientT.trunc_invgamma <- function(y) {
  return(suff.T = cbind(log(y), 1 / y))
}

averageT.trunc_invgamma <- function(y) {
  return(apply(cbind(log(y), 1 / y), 2, mean))
}

#' @export
natural2parameters.trunc_invgamma <- function(eta) {
  # eta: The natural parameters in a inverse gamma distribution
  # returns (shape,rate)
  parms <- c(shape = -eta[1] - 1, rate = -eta[2])
  class(parms) <- class(eta)
  return(parms)
}

#' @export
parameters2natural.trunc_invgamma <- function(parms) {
  # parms: The parameters shape and rate in a beta distribution
  # returns the natural parameters
  eta <- c(shape = -parms[1] - 1, rate = -parms[2])
  class(eta) <- class(parms)
  return(eta)
}

getYseq.trunc_invgamma <- function(y, y.min = 1e-10, y.max = 1, n = 100) {
  # needs chekking
  mean <- mean(y, na.rm = TRUE)
  sd <- var(y, na.rm = TRUE)^0.5
  lo <- max(y.min, mean - 5 * sd, 1e-10)
  hi <- min(y.max, mean + 5 * sd)
  out <- seq(lo, hi, length = n)
  class(out) <- class(y)
  return(out)
}

getGradETinv.trunc_invgamma <- function(eta) {
  # eta: Natural parameter
  # return the inverse of E.T differentiated with respect to eta' : p x p matrix
  A.11 <- sum(1 / (((0:10000) + eta[1] + 1))^2)
  return(A = solve(matrix(c(A.11, -1 / eta[2], -1 / eta[2], (eta[1] + 1) / eta[2]^2, ncol = 2))))
}
