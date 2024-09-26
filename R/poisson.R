## --##--##--##--##--##--##--##--##--##--##--##--##--##
##   Functions related to the Poisson distribution   ##
## --##--##--##--##--##--##--##--##--##--##--##--##--##

#' @param lambda mean and var of "parent" distribution
#' @rdname rtrunc
#' @export
rtruncpois <- rtrunc.poisson <- function(
  n, lambda, a = 0, b = Inf, faster = FALSE
) {
  class(n) <- "trunc_poisson"
  if (faster) {
    family <- gsub("trunc_", "", class(n))
    parms <- mget(ls())[grep("^faster$|^n$|^family$", ls(), invert = TRUE)]
    return(rtrunc_direct(n, family, parms, a, b))
  } else {
    parms <- mget(ls())[grep("^faster$", ls(), invert = TRUE)]
    return(sampleFromTruncated(parms))
  }
}


#' @export
dtrunc.trunc_poisson <- function(y, lambda, eta, a = 0, b = Inf, ...) {
  if (missing(eta)) {
    eta <- parameters2natural.parms_poisson(lambda)
  }
  parm <- exp(eta)
  dens <- rescaledDensities(y, a - 1, b, dpois, ppois, parm)
  return(dens)
}

#' @rdname dtrunc
#' @export
dtruncpois <- dtrunc.trunc_poisson

#' @export
empiricalParameters.trunc_poisson <- function(y, ...) {
  # Returns empirical parameter estimate for lambda
  parms <- c("lambda" = mean(y))
  class(parms) <- "parms_poisson"
  parms
}

#' @method sufficientT trunc_poisson
sufficientT.trunc_poisson <- function(y) {
  y
}

#' @export
natural2parameters.parms_poisson <- function(eta, ...) {
  # eta: The natural parameters in a Poisson distribution
  # returns (mean,sigma)
  if (length(eta) != 1) stop("Eta must be one single number")
  lambda <- c(lambda = exp(eta[[1]]))
  class(lambda) <- class(eta)
  lambda
}

#' @export
parameters2natural.parms_poisson <- function(parms, ...) {
  # parms: The parameter lambda in a Poisson distribution
  # returns the natural parameters
  eta <- prepEta(log(parms), class(parms))
}

#' @method getGradETinv parms_poisson
getGradETinv.parms_poisson <- function(eta, ...) {
  # eta: Natural parameter
  # return the inverse of E.T differentiated with respect to eta
  exp(-eta)
}

#' @method getYseq trunc_poisson
getYseq.trunc_poisson <- function(y, y.min = 0, y.max, n = 100) {
  mean <- mean(y, na.rm = TRUE)
  var.y <- var(y, na.rm = TRUE)
  lo <- max(round(y.min), 0)
  hi <- min(y.max, round(mean + 10 * sqrt(var.y)))
  out <- seq(lo, hi)
  class(out) <- class(y)
  return(out)
}
