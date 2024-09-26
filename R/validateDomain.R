validateDomain <- function(n, ...) {
  UseMethod("validateDomain")
}

#' @method validateDomain trunc_beta
validateDomain.trunc_beta <- function(n, parms, ...) {
  if (parms$shape1 < 0) {
    stop("Invalid parameter domain. shape1 must be non-negative.")
  }
  if (parms$shape2 < 0) {
    stop("Invalid parameter domain. shape2 must be non-negative.")
  }
}

#' @method validateDomain trunc_binomial
validateDomain.trunc_binomial <- function(n, parms, ...) {
  if (parms$size != as.integer(parms$size) || parms$size < 0) {
    stop("Invalid parameter domain. size must be a natural number.")
  }
  if (parms$prob < 0 || parms$prob > 1) {
    stop("Invalid parameter domain. prob must be [0, 1].")
  }
}

#' @method validateDomain trunc_chisq
validateDomain.trunc_chisq <- function(n, parms, ...) {
  if (parms$df != as.integer(parms$df) || parms$df < 0) {
    stop("Invalid parameter domain. df must be a natural positive number.")
  }
}

#' @method validateDomain trunc_contbern
validateDomain.trunc_contbern <- function(n, parms, ...) {
  if (parms$lambda <= 0 || parms$lambda >= 1) {
    stop("Invalid parameter domain. lambda must be (0, 1).")
  }
}

#' @method validateDomain trunc_exp
validateDomain.trunc_exp <- function(n, parms, ...) {
  if (parms$rate <= 0) stop("Invalid parameter domain. rate must be positive.")
}

#' @method validateDomain trunc_gamma
validateDomain.trunc_gamma <- function(n, parms, ...) {
  if (parms$shape <= 0) stop("Invalid parameter domain. shape must be > 0.")
  if (!is.null(parms$rate)) {
      if (parms$rate <= 0) {
        stop("Invalid parameter domain. rate must be > 0.")
      }
  }
  if (!is.null(parms$scale)) {
    if (parms$scale <= 0) {
        stop("Invalid parameter domain. scale must be > 0.")
      }
  }
}

#' @method validateDomain trunc_invgamma
validateDomain.trunc_invgamma <- function(n, parms, ...) {
  if (parms$shape <= 0) stop("Invalid parameter domain. shape must be > 0.")
  if (!is.null(parms$rate)) {
      if (parms$rate <= 0) {
        stop("Invalid parameter domain. rate must be > 0.")
      }
  }
  if (!is.null(parms$scale)) {
    if (parms$scale <= 0) {
        stop("Invalid parameter domain. scale must be > 0.")
      }
  }
}

#' @method validateDomain trunc_invgauss
validateDomain.trunc_invgauss <- function(n, parms, ...) {
  if (parms$m <= 0) stop("Invalid parameter domain. m must be > 0.")
  if (parms$s <= 0) stop("Invalid parameter domain. s must be > 0.")
}

#' @method validateDomain trunc_lognormal
validateDomain.trunc_lognormal <- function(n, parms, ...) {
  if (is.complex(parms$meanlog)) {
    stop("Invalid parameter domain. meanlog must be real.")
  }
  if (parms$sdlog <= 0) stop("Invalid parameter domain. rate must be > 0.")
}

#' @method validateDomain trunc_nbinom
validateDomain.trunc_nbinom <- function(n, parms, ...) {
  if (any(vapply(parms, is.complex, TRUE))) {
    stop("Invalid parameter domain. mean and prob must be real.")
  }
  if (parms$size != as.integer(parms$size) || parms$size < 0) {
    stop("Invalid parameter domain. size must be a natural number.")
  }
  if (!is.na(match("prob", names(parms)))) {
    if (parms$prob < 0 || parms$prob > 1) {
      stop("Invalid parameter domain. prob must be [0, 1].")
    }
  }
}

#' @method validateDomain trunc_normal
validateDomain.trunc_normal <- function(n, parms, ...) {
  if (is.complex(parms$mean)) {
    stop("Invalid parameter domain. mean must be real.")
  }
  if (parms$sd <= 0) stop("Invalid parameter domain. sd must be > 0.")
}

#' @method validateDomain trunc_poisson
validateDomain.trunc_poisson <- function(n, parms, ...) {
  if (parms$lambda < 0) {
    stop("Invalid parameter domain. lambda must be non-negative.")
  }
}
