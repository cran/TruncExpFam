## --##--##--##--##--##--##--##--##--##--##--##--##--##--##--##
##   Functions related to the inverse gamma distribution     ##
## --##--##--##--##--##--##--##--##--##--##--##--##--##--##--##

#' @param shape inverse gamma shape parameter
#' @param rate inverse gamma rate parameter
#' @param scale inverse gamma scale parameter
#' @rdname rtrunc
#' @export
rtruncinvgamma <- rtrunc.invgamma <- function(
  n, shape, rate = 1, scale = 1 / rate, a = 0, b = Inf, faster = FALSE
) {
  class(n) <- "trunc_invgamma"
  if (faster) {
    family <- gsub("trunc_", "", class(n))
    excluded_parms <- "^faster$|^n$|^family$|^rate$|^excluded_parms$"
    parms <- mget(ls())[grep(excluded_parms, ls(), invert = TRUE)]
    return(rtrunc_direct(n, family, parms, a, b))
  } else {
    parms <- mget(ls())[grep("^faster$", ls(), invert = TRUE)]
    return(sampleFromTruncated(parms))
  }
}

#' @export
dtrunc.trunc_invgamma <- function(
  y, shape, rate = 1, scale = 1 / rate, eta, a = 0, b = Inf, ...
) {
  if (missing(eta)) {
    eta <- parameters2natural.parms_invgamma(
      c("shape" = shape, "rate" = rate, "scale" = scale)
    )
  }
  parm <- natural2parameters.parms_invgamma(eta)
  dens <- rescaledDensities(
    y, a, b, dinvgamma, pinvgamma, parm["shape"], parm["rate"]
  )
}

#' @rdname dtrunc
#' @export
dtruncinvgamma <- dtrunc.trunc_invgamma

#' @export
empiricalParameters.trunc_invgamma <- function(y, ...) {
  # Returns  parameter estimates mean and sd
  amean <- mean(y)
  avar <- var(y)
  alpha <- amean^2 / avar + 2
  beta <- (alpha - 1) * amean
  parms <- c(shape = alpha, rate = beta)
  class(parms) <- "parms_invgamma"
  parms
}

#' @method sufficientT trunc_invgamma
sufficientT.trunc_invgamma <- function(y) {
  cbind(log(y), 1 / y)
}

#' @export
natural2parameters.parms_invgamma <- function(eta, ...) {
  # eta: The natural parameters in a inverse gamma distribution
  # returns (shape,rate)
  if (length(eta) != 2) stop("Eta must be a vector of two elements")
  parms <- c("shape" = -eta[[1]] - 1, "rate" = -eta[[2]])
  class(parms) <- class(eta)
  parms
}

#' @export
parameters2natural.parms_invgamma <- function(parms, ...) {
  # parms: The parameters shape and rate in a beta distribution
  # returns the natural parameters
  eta <- c(eta1 = -parms[[1]] - 1, eta2 = -parms[[2]])
  class(eta) <- class(parms)
  eta
}

#' @method getYseq trunc_invgamma
getYseq.trunc_invgamma <- function(y, y.min = 1e-10, y.max = 1, n = 100) {
  # needs chekking
  mean <- mean(y, na.rm = TRUE)
  sd <- var(y, na.rm = TRUE)^0.5
  lo <- max(y.min, mean - 5 * sd, 1e-10)
  hi <- min(y.max, mean + 5 * sd)
  out <- seq(lo, hi, length = n)
  class(out) <- class(y)
  return(out)
}

#' @method getGradETinv parms_invgamma
getGradETinv.parms_invgamma <- function(eta, ...) {
  # eta: Natural parameter
  # return the inverse of E.T differentiated with respect to eta' : p x p matrix
  A.11 <- sum(1 / (0:10000 + eta[1] + 1)^2)
  A.22 <- sum((0:10000 + eta[1] + 1) / eta[2]^2)
  A.12 <- -1 / eta[2]
  inv_A <- matrix(c(A.11, A.12, A.12, A.22), ncol = 2)
  solve(inv_A)
}
