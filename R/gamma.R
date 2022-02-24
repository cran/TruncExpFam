## --##--##--##--##--##--##--##--##--##--##--##--##--##
##   Functions related to the gamma distribution     ##
## --##--##--##--##--##--##--##--##--##--##--##--##--##

#' @param shape shape of "parent" distribution
#' @param rate rate of "parent" distribution
#' @param scale scale of "parent" distribution
#' @rdname rtrunc
#' @export
rtruncgamma <- rtrunc.gamma <- function(n, shape, rate = 1, scale = 1 / rate,
                                        a = 0, b = Inf) {
  class(n) <- "trunc_gamma"
  sampleFromTruncated(mget(ls()))
}

#' @export
dtrunc.trunc_gamma <- function(y, eta, a = 0, b = Inf) {
  parm <- natural2parameters.trunc_gamma(eta)
  dens <- ifelse((y < a) | (y > b), 0, dgamma(y, shape = parm[1], rate = parm[2]))
  if (!missing(a)) {
    F.a <- pgamma(a, shape = parm[1], rate = parm[2])
  } else {
    F.a <- 0
  }
  if (!missing(b)) {
    F.b <- pgamma(b, shape = parm[1], rate = parm[2])
  } else {
    F.b <- 1
  }
  const <- 1 / (F.b - F.a)
  return(dens * const)
}

#' @rdname dtrunc
#' @export
dtruncgamma <- dtrunc.trunc_gamma

#' @export
init.parms.trunc_gamma <- function(y) {
  # Returns  parameter estimates mean and sd
  amean <- mean(y)
  avar <- var(y)
  shp <- amean^2 / avar
  parms <- c(shape = shp, rate = shp / amean)
  class(parms) <- "trunc_gamma"
  return(parms)
}

sufficientT.trunc_gamma <- function(y) {
  return(suff.T = cbind(log(y), y))
}

averageT.trunc_gamma <- function(y) {
  return(apply(cbind(log(y), y), 2, mean))
}

#' @export
natural2parameters.trunc_gamma <- function(eta) {
  # eta: The natural parameters in a gamma distribution
  # returns (shape,rate)
  parms <- c(shape = eta[1] + 1, rate = -eta[2])
  class(parms) <- class(eta)
  return(parms)
}

#' @export
parameters2natural.trunc_gamma <- function(parms) {
  # parms: The parameters shape and rate in a gamma distribution
  # returns the natural parameters
  eta <- c(eta.1 = parms[1] - 1, eta.2 = -parms[2])
  class(eta) <- class(parms)
  return(eta)
}

getYseq.trunc_gamma <- function(y, y.min = 1e-6, y.max, n = 100) {
  # Bør chekkes
  mean <- mean(y, na.rm = TRUE)
  sd <- var(y, na.rm = TRUE)^0.5
  lo <- max(y.min, mean - 5 * sd)
  hi <- min(y.max, mean + 5 * sd)
  out <- seq(lo, hi, length = n)
  class(out) <- class(y)
  return(out)
}

getGradETinv.trunc_gamma <- function(eta) {
  # eta: Natural parameter
  # return the inverse of E.T differentiated with respect to eta' : p x p matrix
  return(A = solve(matrix(c(-1 / eta[1]^2 + dpsi.dx(eta[1]), -1 / eta[2], -1 / eta[2], (eta[1] + 1) / eta[2]^2), ncol = 2)))
}
