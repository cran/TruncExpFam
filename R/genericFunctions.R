#' @title Initialize parameters
#' @description Returns the empirical parameter estimate for a distribution
#' @param y output of rtrunc
#' @examples
#' # Normal distribution
#' sampNorm <- rtrunc(50, mean = 5, sd = 2)
#' init.parms(sampNorm)
#'
#' # Poisson distribution
#' sampPois <- rtrunc(10, lambda = 100, family = "Poisson")
#' init.parms(sampPois)
#' @export
#' @return A vector of parameter estimates for the input sample
init.parms <- function(y) {
  UseMethod("init.parms")
}

#' @title Convert natural parameters to distribution parameters
#' @param eta vector of natural parameters
#' @return A vector of the original distribution parameters
#' @export
#' @examples
#' samp <- rtrunc(n = 100, lambda = 2, family = "Poisson")
#' natural2parameters(init.parms(samp))
natural2parameters <- function(eta) {
  UseMethod("natural2parameters")
}

#' @title Convert distribution parameters to natural parameters
#' @param parms A vector of parameters in a distribution distribution
#' @return A vector containing the natural parameters
#' @export
#' @examples
#' # Poisson distribution
#' samp <- rtrunc(n = 100, lambda = 2, family = "Poisson")
#' parameters2natural(init.parms(samp))
parameters2natural <- function(parms) {
  UseMethod("parameters2natural")
}

sufficientT <- function(y) {
  UseMethod("sufficientT")
}

getYseq <- function(y, y.min, y.max, n) {
  UseMethod("getYseq", y)
}

averageT <- function(y) {
  UseMethod("averageT")
}

getGradETinv <- function(eta, ...) {
  UseMethod("getGradETinv")
}
