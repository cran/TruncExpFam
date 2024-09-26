## --##--##--##--##--##--##--##--##--##--##--##--##--##--##
##   Functions related to the Chi Square distribution    ##
## --##--##--##--##--##--##--##--##--##--##--##--##--##--##

#' @param rate vector of rates
#' @rdname rtrunc
#' @export
rtruncexp <- function(n, rate = 1, a = 0, b = Inf, faster = FALSE) {
  class(n) <- "trunc_exp"
  if (faster) {
    family <- gsub("trunc_", "", class(n))
    parms <- mget(ls())[grep("^faster$|^n$|^family$", ls(), invert = TRUE)]
    return(rtrunc_direct(n, family, parms, a, b))
  } else {
    parms <- mget(ls())[grep("^faster$", ls(), invert = TRUE)]
    return(sampleFromTruncated(parms))
  }
}
rtrunc.exp <- rtruncexp

#' @export
dtrunc.trunc_exp <- function(y, rate = 1, eta, a = 0, b = Inf, ...) {
  if (missing(eta)) {
    eta <- parameters2natural.parms_exp(c("rate" = rate))
  }
  rate <- natural2parameters.parms_exp(eta)
  dens <- rescaledDensities(y, a, b, dexp, pexp, rate)
  return(dens)
}

#' @rdname dtrunc
#' @export
dtruncexp <- dtrunc.trunc_exp

#' @export
empiricalParameters.trunc_exp <- function(y, ...) {
  # Returns empirical parameter estimate for the rate parameter
  parms <- c("rate" = mean(y))
  class(parms) <- "parms_exp"
  parms
}

#' @method sufficientT trunc_exp
sufficientT.trunc_exp <- function(y) {
  y
}

#' @export
natural2parameters.parms_exp <- function(eta, ...) {
  # eta: The natural parameters in an exponential distribution distribution
  # returns rate
  if (length(eta) != 1) stop("Eta must be one single number")
  lambda <- c(rate = -eta[[1]])
  class(lambda) <- class(eta)
  lambda
}

#' @export
parameters2natural.parms_exp <- function(parms, ...) {
  # parms: The parameter lambda in an exponential distribution
  # returns the natural parameters
  eta <- c("eta" = -parms[["rate"]])
  class(eta) <- class(parms)
  eta
}

#' @method getGradETinv parms_exp
getGradETinv.parms_exp <- function(eta, ...) {
  # eta: Natural parameter
  # return the inverse of E.T differentiated with respect to eta
  eta^2
}

#' @method getYseq trunc_exp
getYseq.trunc_exp <- function(y, y.min = 0, y.max, n = 100) {
  mean <- mean(y, na.rm = TRUE)
  var.y <- var(y, na.rm = TRUE)
  lo <- max(round(y.min), 0)
  hi <- min(y.max, round(mean + 10 * sqrt(var.y)))
  out <- seq(lo, hi, length = n)
  class(out) <- class(y)
  return(out)
}
