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
dtrunc.trunc_normal <- function(y, mean = 0, sd = 1, eta, a = -Inf, b = Inf, ...) {
  if (missing(eta)) {
    eta <- parameters2natural.parms_normal(c("mean" = mean, "sd" = sd))
  }
  parm <- natural2parameters.parms_normal(eta)
  dens <- rescaledDensities(y, a, b, dnorm, pnorm, parm[1], parm[2])
  return(dens)
}

#' @rdname dtrunc
#' @export
dtruncnorm <- dtrunc.trunc_normal

#' @export
empiricalParameters.trunc_normal <- function(y, ...) {
  # Returns empirical parameter estimates mean and sd
  parms <- c(mean = mean(y), sd = sqrt(var(y)))
  class(parms) <- "parms_normal"
  return(parms)
}

sufficientT.trunc_normal <- function(y) {
  return(suff.T = cbind(y, y^2))
}

#' @export
natural2parameters.parms_normal <- function(eta, ...) {
  # eta: The natural parameters in a normal distribution
  # returns (mean,sigma)
  if (length(eta) != 2) stop("Eta must be a vector of two elements")
  parms <- c("mean" = -0.5 * eta[[1]] / eta[[2]], "sd" = sqrt(-0.5 / eta[[2]]))
  class(parms) <- class(eta)
  return(parms)
}

#' @export
parameters2natural.parms_normal <- function(parms, ...) {
  # parms: The parameters mean and sd in a normal distribution
  # returns the natural parameters
  eta <- c(eta1 = parms[["mean"]], eta2 = -0.5) / parms[["sd"]]^2
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

getGradETinv.parms_normal <- function(eta, ...) {
  # eta: Natural parameter
  # return the inverse of E.T differentiated with respect to eta' : p x p matrix
  A_inv <- 0.5 * matrix(
    c(
      -1 / eta[2], eta[1] / eta[2]^2,
      eta[1] / eta[2]^2, 1 / eta[2]^2 - eta[1]^2 / eta[2]^3
    ),
    ncol = 2
  )
  return(A = solve(A_inv))
}
