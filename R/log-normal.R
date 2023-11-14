## --##--##--##--##--##--##--##--##--##--##--##--##--##--##
##   Functions related to the Log Normal distribution   ##
## --##--##--##--##--##--##--##--##--##--##--##--##--##--##

#' @param meanlog mean of un-truncated distribution
#' @param sdlog standard deviation of un-truncated distribution
#' @rdname rtrunc
#' @export
rtrunclnorm <- rtrunc.lognormal <- function(n, meanlog, sdlog, a = 0, b = Inf) {
  class(n) <- "trunc_lognormal"
  sampleFromTruncated(mget(ls()))
}

sufficientT.trunc_lognormal <- function(y) {
  return(suff.T = cbind(log(y), log(y)^2))
}

#' @export
dtrunc.trunc_lognormal <- function(
  y, meanlog = 0, sdlog = 1, eta, a = 0, b = Inf, ...
) {
  if (missing(eta)) {
    eta <- parameters2natural.parms_lognormal(c("meanlog" = meanlog, "sdlog" = sdlog))
  }
  parm <- natural2parameters.parms_lognormal(eta)
  dens <- rescaledDensities(y, a, b, dlnorm, plnorm, parm[1], parm[2])
  return(dens)
}

#' @rdname dtrunc
#' @export
dtrunclnorm <- dtrunc.trunc_lognormal

#' @export
empiricalParameters.trunc_lognormal <- function(y, ...) {
  # Y~LN(mean,sigma) => X=log(Y)~N(mean,sigma)
  # Returns empirical parameter estimates for mean and sd
  x <- log(y)
  parms <- c("meanlog" = mean(x), "sdlog" = sqrt(var(x)))
  class(parms) <- "parms_lognormal"
  return(parms)
}

getYseq.trunc_lognormal <- function(y, y.min, y.max, n = 100) {
  x <- log(y)
  mean <- mean(x, na.rm = TRUE)
  sd <- var(x, na.rm = TRUE)^0.5
  lo <- max(y.min, exp(mean - 3.5 * sd))
  hi <- min(y.max, exp(mean + 3.5 * sd))
  out <- seq(lo, hi, length = n)
  class(out) <- class(y)
  return(out)
}

#' @export
natural2parameters.parms_lognormal <- function(eta, ...) {
  if (length(eta) != 2) stop("Eta must be a vector of two elements")
  parms <- c("meanlog" = -0.5 * eta[[1]] / eta[[2]], "sdlog" = sqrt(-0.5 / eta[[2]]))
  class(parms) <- class(eta)
  return(parms)
}

#' @export
parameters2natural.parms_lognormal <- function(parms, ...) {
  eta <- c(eta1 = parms[["meanlog"]], eta2 = -0.5) / parms[["sdlog"]]^2
  class(eta) <- class(parms)
  return(eta)
}

getGradETinv.parms_lognormal <- function(eta, ...) {
  getGradETinv.parms_normal(eta, ...)
}
