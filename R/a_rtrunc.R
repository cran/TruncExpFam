#' @title The Truncated Exponential Family
#' @description Random generation for the truncated exponential family distributions. Please ferer to the "Details" and "Examples" section for more information on how to use this function.
#' @param n sample size
#' @param a point of left truncation
#' @param b point of right truncation
#' @param family distribution family to use
#' @param ... individual arguments to each distribution
#' @return A sample of size n drawn from a truncated distribution
#' @note The current sample-generating algorithm may be slow if the distribution
#' is largely represented by low-probability values. This will be fixed soon.
#' Please follow \url{https://github.com/ocbe-uio/TruncExpFam/issues/72} for
#' details.
#' @author René Holst, Waldir Leôncio
#' @details One way to use this function is by calling the `rtrunc`
#' generic with the `family` parameter of your choice. You can also
#' specifically call one of the methods (e.g. `rtrunc.poisson(10, lambda=3)`
#' instead of `rtrunc(10, family="poisson", lambda=3)). The latter is more
#' flexible (i.e., easily programmable) and more robust (i.e., it contains
#' better error handling and validation procedures), while the former better
#' conforms with the nomenclature from other distribution-related functions in
#' the \code{stats} package.
#' @importFrom methods new
#' @return vector of one of the \code{rtrunc_*} classes containing the sample
#' elements, as well as some attributes related to the chosen distribution.
#' @examples
#' # Truncated binomial distribution
#' sample.binom <- rtrunc(100, family = "binomial", prob = 0.6, size = 20, a = 4, b = 10)
#' sample.binom
#' plot(table(sample.binom), ylab = "Frequency", main = "Freq. of sampled values")
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
rtrunc <- function(n, family = "gaussian", ...) {

  # Validating ---------------------------------------------------------------
  family <- tolower(family)
  validateFamilyName(family)

  # Determining object class -------------------------------------------------
  parms <- list(...)
  trunc_class <- genrtruncClass(n, family, names(parms))
  extra_n <- 1 # to generate extra observations to complete n from input
  class(extra_n) <- class(n) <- trunc_class

  # Generating sample --------------------------------------------------------
  sample <- rtrunc.generic(n, ...)
  saved_attributes <- attributes(sample)
  while (length(sample) != n) {
    new_obs <- rtrunc.generic(extra_n, ...)
    sample <- c(sample, new_obs)
    class(sample) <- class(new_obs)
  }

  # Attaching attributes -----------------------------------------------------
  attributes(sample) <- saved_attributes

  # Returning sampled elements -----------------------------------------------
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

validateFamilyName <- function(family) {
  family <- tolower(family)
  if (!(family %in% valid_distros)) {
    stop(
      "Invalid distribution family. Please choose from the list below:\n",
      paste(valid_distros, collapse = ", ")
    )
  }
}

#' @title Validate family parameters
#' @description Checks if a combination of distribution family and parameters is
#' valid.
#' @param family character with family distribution name
#' @param parms character vector with distribution parameter names
#' @param verbose print intermediate messages?
#' @return list telling if family-parm combo is valid + the family name
#' @author Waldir Leoncio
validateFamilyParms <- function(family, parms, verbose = FALSE) {
  matched <- list(family = FALSE, parameters = FALSE)
  families <- grep(family, names(valid_fam_parm))
  for (fam in families) {
    if (any(family == valid_fam_parm[[fam]]$family)) {
      if (verbose) message("Matched family: ", family)
      matched$family <- TRUE
      family <- valid_fam_parm[[fam]]$family[1] # use standard family name
      parms_text <- paste(parms, collapse = ", ")
      parms_expected <- valid_fam_parm[[fam]]$parms
      if (all(sort(parms) == sort(parms_expected))) {
        if (verbose) message("Matched parameters: ", parms_text)
        matched$parameters <- TRUE
      }
    }
  }
  if ("parms_expected" %in% ls() & !matched$parameters) {
    parms_expected_text <- paste(unlist(parms_expected), collapse = ", ")
    stop(
      "The {", parms_text, "} ",
      "parameter set does not match the ", family, " family. ",
      "Expected set of parameters: {", parms_expected_text, "}. ",
      "Please change the family to match the expected ",
      "parameters or use a different family."
    )
  }
  return(list(is_valid = all(unlist(matched)), family_name = family))
}
