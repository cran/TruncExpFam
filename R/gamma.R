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
  if (!missing(rate) && !missing(scale)) {
    stop("specify 'rate' or 'scale' but not both")
  }
  class(n) <- "trunc_gamma"
  sampleFromTruncated(mget(ls()))
}

#' @export
dtrunc.trunc_gamma <- function(
  y, shape, rate = 1, scale = 1 / rate, eta, a = 0, b = Inf, ...
) {
  if (missing(eta)) {
    eta <- parameters2natural.parms_gamma(c("shape" = shape, "rate" = rate, "scale" = scale))
  }
  parm <- natural2parameters.parms_gamma(eta)
  dens <- rescaledDensities(y, a, b, dgamma, pgamma, parm["shape"], parm["rate"])
  return(dens)
}

#' @rdname dtrunc
#' @export
dtruncgamma <- dtrunc.trunc_gamma

#' @export
empiricalParameters.trunc_gamma <- function(y, ...) {
  # Returns  parameter estimates mean and sd
  amean <- mean(y)
  avar <- var(y)
  shp <- amean^2 / avar
  parms <- c(shape = shp, rate = shp / amean)
  class(parms) <- "parms_gamma"
  return(parms)
}

sufficientT.trunc_gamma <- function(y) {
  return(suff.T = cbind(log(y), y))
}

#' @export
natural2parameters.parms_gamma <- function(eta, ...) {
  # eta: The natural parameters in a gamma distribution
  # returns (shape,rate)
  if (length(eta) != 2) stop("Eta must be a vector of two elements")
  parms <- c("shape" = eta[[1]] + 1, "rate" = -eta[[2]])
  class(parms) <- class(eta)
  return(parms)
}

#' @export
parameters2natural.parms_gamma <- function(parms, ...) {
  # parms: The parameters shape and rate in a gamma distribution
  # returns the natural parameters
  if (all(c("shape", "rate") %in% names(parms))) {
    eta <- c(eta1 = parms[["shape"]] - 1, eta2 = -parms[["rate"]])
  } else {
    eta <- c(eta1 = parms[["shape"]] - 1, eta2 = -1 / parms[["scale"]])
  }
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
  out <- out[out > 0] # prevents NaN as sufficient statistics
  class(out) <- class(y)
  return(out)
}

getGradETinv.parms_gamma <- function(eta, ...) {
  # eta: Natural parameter
  # return the inverse of E.T differentiated with respect to eta' : p x p matrix
  dpsi.dx <- function(x, k = 10000) {
    # Returns the derivative of the psi function (removed)
    sum((1 / ((0:k) + x))^2)
  }
  A_inv <- matrix(
    c(
      -1 / eta[1]^2 + dpsi.dx(eta[1]), -1 / eta[2],
      -1 / eta[2], (eta[1] + 1) / eta[2]^2
    ),
    ncol = 2
  )
  return(A = solve(A_inv))
}
