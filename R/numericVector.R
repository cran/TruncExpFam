#' @title Extract parameters
#' @param y Numeric vector containing observations from a random variable
#' @param family Distribution family to assume for `y`
#' @param natural Should output be in terms of the natural parameter eta?
#' @param ... arguments passed to [empiricalParameters()]
#' @export
#' @examples
#' # Some random data
#' x <- c(
#'   4, 3, 6, 3, 3, 3, 3, 4, 3, 2, 3, 0, 4, 2, 0, 1, 4, 3, 0, 0, 2, 3, 0, 3, 7,
#'   2, 1, 1, 2, 3, 2, 3, 3, 3, 2, 2, 2, 0, 2, 0, 2, 1, 0, 2, 3, 1, 0, 4, 2, 2,
#'   0, 1, 1, 1, 2, 2, 3, 1, 3, 1, 1, 0, 3, 3, 2, 0, 2, 2, 3, 0, 2, 1, 0, 0, 1,
#'   0, 2, 4, 2, 3, 3, 0, 1, 0, 5, 2, 4, 2, 7, 4, 4, 1, 2, 4, 3, 2, 4, 3, 1, 3
#' )
#'
#' # Extracting parameters under different distribution assumptions
#' empiricalParameters(x, family = "normal")
#' empiricalParameters(x, family = "normal", natural = TRUE)
#' empiricalParameters(x, family = "binomial", nsize = max(x))
#' empiricalParameters(x, family = "poisson", natural = FALSE)
#' empiricalParameters(x, family = "poisson", natural = TRUE)
empiricalParameters.numeric <- function(
  y, family = "gaussian", natural = FALSE, ...
) {

  # Assigning trunc family, if necessary
  class(y) <- paste0("trunc_", useStandardFamilyName(family))

  # Estimating and validating parameters
  parms <- empiricalParameters(y, ...)
  validateSupport(y, as.list(parms), ...)
  if (natural) parms <- parameters2natural(parms)
  return(parms)
}

#' @export
parameters2natural.numeric <- function(parms, family = "gaussian", ...) {
  # Validation and reformatting
  family <- useStandardFamilyName(family)
  if (is.null(names(parms))) {
    stop("Please provide a named parameter vector")
  } else {
    validateFamilyParms(family, names(parms))
  }
  class(parms) <- paste0("parms_", family)

  # Dispatching
  return(parameters2natural(parms))
}

#' @export
natural2parameters.numeric <- function(eta, family = "gaussian", ...) {
  # Validation and reformatting
  family <- useStandardFamilyName(family)
  validateNaturalParms(names(eta))
  class(eta) <- paste0("parms_", family)

  # Dispatching
  return(natural2parameters(eta))
}
