## --##--##--##--##--##--##--##--##--##--##--##--##--##--##
##   Functions related to the Chi Square distribution    ##
## --##--##--##--##--##--##--##--##--##--##--##--##--##--##

#' @param df degrees of freedom for "parent" distribution
#' @rdname rtrunc
#' @export
rtruncchisq <- rtrunc.chisq <- function(n, df, a = 0, b = Inf, faster = FALSE) {
  class(n) <- "trunc_chisq"
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
dtrunc.trunc_chisq <- function(y, df, eta, a = 0, b = Inf, ...) {
  if (missing(eta)) {
    eta <- parameters2natural.parms_chisq(c("df" = df))
  }
  df <- natural2parameters.parms_chisq(eta)
  dens <- rescaledDensities(y, a, b, dchisq, pchisq, df)
  return(dens)
}

#' @rdname dtrunc
#' @export
dtruncchisq <- dtrunc.trunc_chisq

#' @export
empiricalParameters.trunc_chisq <- function(y, ...) {
  # Returns empirical parameter estimate for df
  parms <- c("df" = mean(y))
  class(parms) <- "parms_chisq"
  parms
}

#' @method sufficientT trunc_chisq
sufficientT.trunc_chisq <- function(y) {
  log(y)
}

#' @export
natural2parameters.parms_chisq <- function(eta, ...) {
  # eta: The natural parameters in a Chi Square distribution
  # returns df
  if (length(eta) != 1) stop("Eta must be one single number")
  df <- c(df = 2 * (eta[[1]] + 1))
  class(df) <- class(eta)
  df
}

#' @export
parameters2natural.parms_chisq <- function(parms, ...) {
  # parms: The parameter lambda in a Chi Square distribution
  # returns the natural parameters
  eta <- prepEta(parms / 2 - 1, class(parms))
}

#' @method getGradETinv parms_chisq
getGradETinv.parms_chisq <- function(eta, ...) {
  # eta: Natural parameter
  # return the inverse of E.T differentiated with respect to eta
  1 / sum(1 / (as.vector(eta) + (1:1e6))^2)
}

#' @method getYseq trunc_chisq
getYseq.trunc_chisq <- function(y, y.min = 0, y.max, n = 100) {
  mean <- mean(y, na.rm = TRUE)
  var.y <- var(y, na.rm = TRUE)
  lo <- max(round(y.min), 0)
  hi <- min(y.max, round(mean + 10 * sqrt(var.y)))
  out <- seq(lo, hi, length = n)
  out <- out[out > 0] # prevents NaN as sufficient statistics
  class(out) <- class(y)
  return(out)
}
