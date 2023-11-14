## --##--##--##--##--##--##--##--##--##--##--##--##--##--##--##
##   Functions related to the inverse gamma distribution     ##
## --##--##--##--##--##--##--##--##--##--##--##--##--##--##--##

#' @param shape inverse gamma shape parameter
#' @param rate inverse gamma rate parameter
#' @param scale inverse gamma scale parameter
#' @rdname rtrunc
#' @export
rtruncinvgamma <- rtrunc.invgamma <- function(
  n, shape, rate = 1, scale = 1 / rate, a = 0, b = Inf
) {
  class(n) <- "trunc_invgamma"
  sampleFromTruncated(mget(ls()))
}

#' @export
dtrunc.trunc_invgamma <- function(
  y, shape, rate = 1, scale = 1 / rate, eta, a = 0, b = Inf, ...
) {
  if (missing(eta)) {
    eta <- parameters2natural.parms_invgamma(c("shape" = shape, "rate" = rate, "scale" = scale))
  }
  parm <- natural2parameters.parms_invgamma(eta)
  dens <- rescaledDensities(y, a, b, dinvgamma, pinvgamma, parm["shape"], parm["rate"])
}

#' @rdname dtrunc
#' @export
dtruncinvgamma <- dtrunc.trunc_invgamma

#' @export
empiricalParameters.trunc_invgamma <- function(y, ...) {
  # Returns  parameter estimates mean and sd
  amean <- mean(y)
  avar <- var(y)
  alpha <- amean^2 / avar + 2
  beta <- (alpha - 1) * amean
  parms <- c(shape = alpha, rate = beta)
  class(parms) <- "parms_invgamma"
  return(parms)
}

sufficientT.trunc_invgamma <- function(y) {
  return(suff.T = cbind(log(y), 1 / y))
}

#' @export
natural2parameters.parms_invgamma <- function(eta, ...) {
  # eta: The natural parameters in a inverse gamma distribution
  # returns (shape,rate)
  if (length(eta) != 2) stop("Eta must be a vector of two elements")
  parms <- c("shape" = -eta[[1]] - 1, "rate" = -eta[[2]])
  class(parms) <- class(eta)
  return(parms)
}

#' @export
parameters2natural.parms_invgamma <- function(parms, ...) {
  # parms: The parameters shape and rate in a beta distribution
  # returns the natural parameters
  eta <- c(eta1 = -parms[[1]] - 1, eta2 = -parms[[2]])
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

getGradETinv.parms_invgamma <- function(eta, ...) {
  # eta: Natural parameter
  # return the inverse of E.T differentiated with respect to eta' : p x p matrix
  A.11 <- sum(1 / (0:10000 + eta[1] + 1)^2)
  A.22 <- sum((0:10000 + eta[1] + 1) / eta[2]^2)
  A.12 <- -1 / eta[2]
  inv_A <- matrix(c(A.11, A.12, A.12, A.22), ncol = 2)
  return(A = solve(inv_A))
}
