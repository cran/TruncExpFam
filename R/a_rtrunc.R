#' @title The Truncated Exponential Family
#' @description Random generation for the truncated exponential family
#' distributions. Please ferer to the "Details" and "Examples" section for more
#' information on how to use this function.
#' @param n sample size
#' @param a point of left truncation
#' @param b point of right truncation
#' @param family distribution family to use
#' @param faster if \code{TRUE}, samples directly from the truncated
#' distribution (more info in details)
#' @param ... individual arguments to each distribution
#' @return A sample of size n drawn from a truncated distribution
#' @note The current sample-generating algorithm may be slow if the distribution
#' is largely represented by low-probability values. This will be fixed soon.
#' Please follow \url{https://github.com/ocbe-uio/TruncExpFam/issues/72} for
#' details.
#' @author René Holst, Waldir Leôncio
#' @details One way to use this function is by calling the `rtrunc`
#' generic with the `family` parameter of your choice. You can also
#' specifically call one of the methods (e.g. `rtruncpois(10, lambda=3)`
#' instead of `rtrunc(10, family="poisson", lambda=3)). The latter is more
#' flexible (i.e., easily programmable) and more robust (i.e., it contains
#' better error handling and validation procedures), while the former better
#' conforms with the nomenclature from other distribution-related functions in
#' the \code{stats} package.
#'
#' Setting \code{faster=TRUE} uses a new algorithm that samples directly from
#' the truncated distribution, as opposed to the old algorithm that samples
#' from the untruncated distribution and then truncates the result. The
#' advantage of the new algorithm is that it is way faster than the old one,
#' particularly for highly-truncated distributions. On the other hand, the
#' sample for untruncated distributions called through `rtrunc()` will no longer
#' match their \code{stats}-package counterparts for the same seed.
#'
#' @return vector of one of the \code{rtrunc_*} classes containing the sample
#' elements, as well as some attributes related to the chosen distribution.
#' @examples
#' # Truncated binomial distribution
#' sample.binom <- rtrunc(
#'   100, family = "binomial", prob = 0.6, size = 20, a = 4, b = 10
#' )
#' sample.binom
#' plot(
#'   table(sample.binom), ylab = "Frequency", main = "Freq. of sampled values"
#' )
#'
#' # Truncated Log-Normal distribution
#' sample.lognorm <- rtrunc(
#'   n = 100, family = "lognormal", meanlog = 2.5, sdlog = 0.5, a = 7
#' )
#' summary(sample.lognorm)
#'
#' hist(
#'   sample.lognorm,
#'   nclass = 35, xlim = c(0, 60), freq = FALSE,
#'   ylim = c(0, 0.15)
#' )
#'
#' # Normal distribution
#' sample.norm <- rtrunc(n = 100, mean = 2, sd = 1.5, a = -1)
#' head(sample.norm)
#' hist(sample.norm, nclass = 25)
#'
#' # Gamma distribution
#' sample.gamma <- rtrunc(n = 100, family = "gamma", shape = 6, rate = 2, a = 2)
#' hist(sample.gamma, nclass = 15)
#'
#' # Poisson distribution
#' sample.pois <- rtrunc(n = 10, family = "poisson", lambda = 10, a = 4)
#' sample.pois
#' plot(table(sample.pois))
#' @export
rtrunc <- function(n, family = "gaussian", faster = FALSE, ...) {
  # This is **not a generic function**, so it doesn't do any class-dispatching.
  # Instead, it relies on the "family" argument for dispatching. Internally,
  # however, this function does use rtrunc.generic() for S3-dispatching
  # according to the class of n.

  # Validating ---------------------------------------------------------------
  family <- tolower(family)
  validateFamilyName(family)

  # Determining object class -------------------------------------------------
  parms <- list(...)
  trunc_class <- genrtruncClass(n, family, names(parms))
  extra_n <- 1 # to generate extra observations to complete n from input
  class(extra_n) <- class(n) <- trunc_class

  # Generating sample --------------------------------------------------------
  if (faster) {
    sample <- rtrunc_direct(n, family, ...)
  } else {
    sample <- rtrunc.generic(n, ...)
  }
  return(sample)
}

rtrunc.generic <- function(n, ...) {
  UseMethod("rtrunc", n)
}

#' @title Generates an rtrunc-dispatchable class
#' @description Matches a list of arguments to an rtrunc method
#' @param n sample size
#' @param family distribution family
#' @param parms list of parameters passed to rtrunc (through the `...` element)
#' @return A character string.
#' @author Waldir Leoncio
genrtruncClass <- function(n, family, parms) {

  # Dropping a and b (parameters not used for validating) -- #
  parms <- parms[!(parms %in% c("a", "b"))]

  # Validating --------------------------------------------- #
  validation_family_parms <- validateFamilyParms(family, parms)
  if (validation_family_parms$is_valid) {
    family <- validation_family_parms$family_name
    return(family)
  }
}
