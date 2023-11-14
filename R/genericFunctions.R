#' @title Calculate empirical parameters
#' @description Returns the empirical parameter estimate for a distribution
#' @param y output of rtrunc
#' @param ... other arguments passed to methods
#' @return A vector of parameter estimates for the input sample
#' @export
#' @examples
#' # Normal distribution
#' sampNorm <- rtrunc(50, mean = 5, sd = 2)
#' empiricalParameters(sampNorm)
#'
#' # Poisson distribution
#' sampPois <- rtrunc(10, lambda = 100, family = "Poisson")
#' empiricalParameters(sampPois)
empiricalParameters <- function(y, ...) {
  UseMethod("empiricalParameters")
}

#' @title Convert natural parameters to distribution parameters
#' @param eta vector of natural parameters
#' @param ... other arguments passed to methods
#' @return A vector of the original distribution parameters
#' @export
#' @seealso [parameters2natural()]
#' @examples
#' samp <- rtrunc(n = 100, lambda = 2, family = "Poisson")
#' lambda_hat <- empiricalParameters(samp)
#' eta_hat <- parameters2natural(lambda_hat)
#' natural2parameters(eta_hat)  # yields back lambda
natural2parameters <- function(eta, ...) {
  UseMethod("natural2parameters")
}

#' @title Convert distribution parameters to natural parameters
#' @param parms A vector of parameters in a distribution distribution
#' @return A vector containing the natural parameters
#' @param ... other arguments passed to methods
#' @export
#' @seealso [natural2parameters()]
#' @examples
#' # Poisson distribution
#' samp <- rtrunc(n = 100, lambda = 2, family = "Poisson")
#' parameters2natural(empiricalParameters(samp))
parameters2natural <- function(parms, ...) {
  UseMethod("parameters2natural")
}

sufficientT <- function(y) {
  UseMethod("sufficientT")
}

getYseq <- function(y, y.min, y.max, n) {
  UseMethod("getYseq", y)
}

getGradETinv <- function(eta, ...) {
  UseMethod("getGradETinv")
}

#' @title Probability Density Function
#' @description Calculates the PDF for a given truncated distribution
#' @note Either the common or the natural parameters must be provided.
#' @param y output from rtrunc or any valid numeric value(s).
#' @export
#' @return The density of y for the given values of the \code{eta} parameter.
#' @examples
#' # Using the output of rtrunc
#' y <- rtrunc(50, mean = 5, sd = 2)
#' dtrunc(y, eta = c(0, -1))
#'
#' # Directly-inputting values
#' dtruncnorm(y = c(5, 0, -10), eta = c(0, -0.05))
dtrunc <- function(y, ...) {
  UseMethod("dtrunc", y)
}
