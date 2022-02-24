## --##--##--##--##--##--##--##--##--##--##--##--##--##--##
##   Functions related to the Chi Square distribution    ##
## --##--##--##--##--##--##--##--##--##--##--##--##--##--##

#' @importFrom stats dchisq pchisq rchisq
#' @param df degrees of freedom for "parent" distribution
#' @rdname rtrunc
#' @export
rtruncchisq <- rtrunc.chisq <- function(n, df, a = 0, b = Inf) {
  class(n) <- "trunc_chisq"
  sampleFromTruncated(mget(ls()))
}

#' @export
dtrunc.trunc_chisq <- function(y, eta, a = 0, b = Inf) {
  df <- natural2parameters.trunc_chisq(eta)
  dens <- ifelse((y <= a) | (y > b), 0, dchisq(y, df = df))
  if (!missing(a)) {
    F.a <- pchisq(a, df)
  } else {
    F.a <- 0
  }
  if (!missing(b)) {
    F.b <- pchisq(b, df)
  } else {
    F.b <- 1
  }
  return(dens / (F.b - F.a))
}

#' @rdname dtrunc
#' @export
dtruncchisq <- dtrunc.trunc_chisq

#' @export
init.parms.trunc_chisq <- function(y) {
  # Returns empirical parameter estimate for df
  parms <- mean(y)
  class(parms) <- "trunc_chisq"
  return(parms)
}

sufficientT.trunc_chisq <- function(y) {
  return(suff.T = log(y))
}

averageT.trunc_chisq <- function(y) {
  return(mean(y))
}

#' @export
natural2parameters.trunc_chisq <- function(eta) {
  # eta: The natural parameters in a Chi Square distribution
  # returns df
  df <- c(parms = 2 * (eta + 1))
  class(df) <- class(eta)
  return(df)
}

#' @export
parameters2natural.trunc_chisq <- function(parms) {
  # parms: The parameter lambda in a Chi Square distribution
  # returns the natural parameters
  eta <- parms / 2 - 1
  class(eta) <- class(parms)
  return(eta)
}

getGradETinv.trunc_chisq <- function(eta) {
  # eta: Natural parameter
  # return the inverse of E.T differentiated with respect to eta
  return(A = 1 / sum(1 / (eta + (1:1000000))^2))
}

getYseq.trunc_chisq <- function(y, y.min = 0, y.max, n = 100) {
  mean <- mean(y, na.rm = TRUE)
  var.y <- var(y, na.rm = TRUE)
  lo <- max(round(y.min), 0)
  hi <- min(y.max, round(mean + 10 * sqrt(var.y)))
  out <- seq(lo, hi, length = n)
  class(out) <- class(y)
  return(out)
}
