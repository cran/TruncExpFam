## --##--##--##--##--##--##--##--##--##--##--##--##--##
##   Functions related to the Binomial distribution  ##
## --##--##--##--##--##--##--##--##--##--##--##--##--##

#' @param size number of trials
#' @param prob probability of success on each trial
#' @rdname rtrunc
#' @export
rtruncbinom <- function(n, size, prob, a = 0, b = size, faster = FALSE) {
  class(n) <- "trunc_binomial"
  if (faster) {
    family <- gsub("trunc_", "", class(n))
    parms <- mget(ls())[grep("^faster$|^n$|^family$", ls(), invert = TRUE)]
    return(rtrunc_direct(n, family, parms, a, b))
  } else {
    parms <- mget(ls())[grep("^faster$", ls(), invert = TRUE)]
    return(sampleFromTruncated(parms))
  }
}
rtrunc.binomial <- rtruncbinom

#' @export
dtrunc.trunc_binomial <- function(
  y, size, prob, eta, a = 0, b = attr(y, "parameters")$size, ...
) {
  if (missing(eta)) {
    eta <- parameters2natural.parms_binomial(c("size" = size, "prob" = prob))
  }
  nsize <- attr(y, "parameters")$size
  proba <- 1 / (1 + exp(-eta))
  dens <- ifelse((y < a) | (y > b), 0, dbinom(y, size = nsize, prob = proba))
  F.a <- pbinom(a - 1L, size = nsize, prob = proba) # -1 bc a = 0 is no trunc
  F.b <- pbinom(b, size = nsize, prob = proba)
  dens <- dens / (F.b - F.a)
  attributes(dens) <- attributes(y)
  return(dens)
}

#' @rdname dtrunc
#' @param ... size
#' @export
dtruncbinom <- dtrunc.trunc_binomial

#' @export
empiricalParameters.trunc_binomial <- function(y, size, ...) {
  # Returns empirical parameter estimates
  if (missing(size)) {
    size <- max(y)
  }
  parms <- c("size" = size, "prob" = mean(y) / size)
  class(parms) <- "parms_binomial"
  parms
}

#' @method sufficientT trunc_binomial
sufficientT.trunc_binomial <- function(y) {
  suff.T <- y
}

#' @export
natural2parameters.parms_binomial <- function(eta, ...) {
  # eta: The natural parameters in a binomial distribution
  # returns (p)
  if (length(eta) != 1) stop("Eta must be one single number")
  p <- c(prob = 1 / (1 + exp(-eta[[1]])))
  class(p) <- class(eta)
  p
}

#' @export
parameters2natural.parms_binomial <- function(parms, ...) {
  # parms: The probability parameter p in a binomial distribution
  # returns the natural parameters
  prob <- parms[["prob"]]
  eta <- prepEta(log(prob / (1 - prob)), class(parms))
  attr(eta, "nsize") <- parms[["size"]]
  eta
}

#' @method getGradETinv parms_binomial
getGradETinv.parms_binomial <- function(eta, ...) {
  # eta: Natural parameter
  # return the inverse of E.T differentiated with respect to eta
  nsize <- attr(eta, "nsize")
  exp.eta <- exp(eta)
  ((1 + exp.eta)^2 / exp.eta) / nsize
}

#' @method getYseq trunc_binomial
getYseq.trunc_binomial <- function(y, y.min = 0, y.max, n = 100) {
  nsize <- attr(y, "parameters")$size
  y.lo <- round(y.min)
  y.hi <- round(y.max)
  lo <- max(y.lo, 0)
  hi <- min(y.hi, nsize)
  out <- seq(lo, hi)
  attributes(out) <- attributes(y)
  return(out)
}
