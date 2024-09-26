## --##--##--##--##--##--##--##--##--##--##--##--##--##--##--##--##
##   Functions related to the Negative Binomial distribution    ##
## --##--##--##--##--##--##--##--##--##--##--##--##--##--##--##--##

#' @param size target for number of successful trials,
#' or dispersion parameter (the shape parameter of the gamma mixing
#' distribution). Must be strictly positive, need not be integer.
#' @param prob probability of success on each trial
#' @param mu alternative parametrization via mean
#' @rdname rtrunc
#' @export
rtruncnbinom <- function(n, size, prob, mu, a = 0, b = Inf, faster = FALSE) {
  class(n) <- "trunc_nbinom"
  if (missing(prob)) {
    prob <- (size) / (size + mu)
    mu <- ""
  }
  if (faster) {
    family <- gsub("trunc_", "", class(n))
    parms <- mget(ls())[grep("^faster$|^n$|^family$|^mu$", ls(), invert = TRUE)]
    return(rtrunc_direct(n, family, parms, a, b))
  } else {
    parms <- mget(ls())[grep("^faster$", ls(), invert = TRUE)]
    return(sampleFromTruncated(parms))
  }
}
rtrunc.nbinom <- rtruncnbinom

#' @rdname dtrunc
#' @param ... size
#' @export
dtrunc.trunc_nbinom <- function(
  y, size, prob, eta, a = 0, b = Inf, ...
) {
  if (missing(eta)) {
    eta <- parameters2natural.parms_nbinom(c("size" = size, "prob" = prob))
  }
  nsize <- attr(y, "parameters")$size
  proba <- attr(y, "parameters")$prob
  dens <- rescaledDensities(y, a, b, dnbinom, pnbinom, nsize, proba)
  return(dens)
}

#' @rdname dtrunc
#' @param ... size
#' @export
dtruncnbinom <- dtrunc.trunc_nbinom

#' @export
#' @rdname dtrunc
dtruncnbinom <- dtrunc.trunc_nbinom

#' @export
empiricalParameters.trunc_nbinom <- function(y, r, k, ...) {
  # Returns empirical parameter estimate for lambda
  if (missing(r) || missing(k)) {
    parms <- c("mean" = mean(y))
  } else {
    parms <- c("size" = r, "prob" = (r - 1) / (r + k - 1))
  }
  class(parms) <- "parms_nbinom"
  parms
}

#' @method sufficientT trunc_nbinom
sufficientT.trunc_nbinom <- function(y) {
  suff.T <- y
}

#' @export
natural2parameters.parms_nbinom <- function(eta, ...) {
  # eta: The natural parameters in a negative binomial distribution
  p <- c(mean = exp(eta))
  class(p) <- class(eta)
  p
}

#' @export
parameters2natural.parms_nbinom <- function(parms, ...) {
  # parms: The p parameter in a negative binomial distribution
  # returns the natural parameters
  if (all(names(parms) == c("size", "prob"))) {
    mean <- parms[["size"]] * (1 - parms[["prob"]]) / parms[["prob"]]
  } else {
    mean <- parms[["mean"]]
  }
  eta <- prepEta(log(mean), class(parms))
}

#' @method getGradETinv parms_nbinom
getGradETinv.parms_nbinom <- function(eta, r = 1e3, ...) {
  # eta: Natural parameter
  # return the inverse of E.T differentiated with respect to eta
  p <- exp(eta)
  r <- exp(r)
  A <- (1 - p) ^ 2 / (r * p)
  return(A)
}

#' @method getYseq trunc_nbinom
getYseq.trunc_nbinom <- function(y, y.min = 0, y.max, n = 100) {
  mean <- mean(y, na.rm = TRUE)
  var.y <- var(y, na.rm = TRUE)
  lo <- max(round(y.min), 0)
  hi <- min(y.max, round(mean + 10 * sqrt(var.y)))
  out <- seq(lo, hi)
  attributes(out) <- attributes(y)
  return(out)
}
