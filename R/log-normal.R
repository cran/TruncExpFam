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

averageT.trunc_lognormal <- function(y) {
  return(apply(sufficientT.trunc_lognormal(y), 2, mean))
}

#' @export
dtrunc.trunc_lognormal <- function(y, eta, a = 0, b = Inf) {
  parm <- natural2parameters.trunc_normal(eta)
  dens <- ifelse((y < a) | (y > b), 0, dlnorm(y, meanlog = parm[1], sdlog = parm[2]))
  if (!missing(a)) {
    F.a <- plnorm(a, parm[1], parm[2])
  } else {
    F.a <- 0
  }
  if (!missing(b)) {
    F.b <- plnorm(b, parm[1], parm[2])
  } else {
    F.b <- 1
  }
  const <- 1 / (F.b - F.a)
  return(dens  * const)
}

#' @rdname dtrunc
#' @export
dtrunclnorm <- dtrunc.trunc_lognormal

#' @export
init.parms.trunc_lognormal <- function(y) {
  # Y~LN(mean,sigma) => X=log(Y)~N(mean,sigma)
  # Returns empirical parameter estimates for mean and sd
  x <- log(y)
  parms <- c(mean = mean(x), sd = sqrt(var(x)))
  class(parms) <- "trunc_lognormal"
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
natural2parameters.trunc_lognormal <- function(eta) {
  natural2parameters.trunc_normal(eta)
}

#' @export
parameters2natural.trunc_lognormal <- function(parms) {
  parameters2natural.trunc_normal(parms)
}

getGradETinv.trunc_lognormal <- function(eta) {
  getGradETinv.trunc_normal(eta)
}
