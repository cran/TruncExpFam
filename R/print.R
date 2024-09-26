#' @title Print sample from truncated distribution
#' @description Special printing methods for trunc_* classes.
#' @param x object to print
#' @param details if `FALSE` (default), hides the attributes of `x`
#' @param ... other arguments passed to [base::print.default()]
#' @return `x` with or without its attributes
#' @author Waldir Leoncio
#' @aliases print.trunc_beta print.trunc_binomial print.trunc_chisq
#' print.trunc_contbern print.trunc_exp print.trunc_gamma print.trunc_invgamma
#' print.trunc_invgauss print.trunc_lognormal print.trunc_nbinom
#' print.trunc_poisson print.trunc_normal print.trunc
#' @export
print.trunc <- function(x, details = FALSE, ...) {
  if (!details) {
    attributes(x)[(names(attributes(x)) != "names")] <- NULL
  }
  print.default(x, ...)
}

#' @export
print.trunc_beta <- function(x, details = FALSE, ...) {
  print.trunc(x, details, ...)
}

#' @export
print.trunc_binomial <- function(x, details = FALSE, ...) {
  print.trunc(x, details, ...)
}

#' @export
print.trunc_chisq <- function(x, details = FALSE, ...) {
  print.trunc(x, details, ...)
}

#' @export
print.trunc_contbern <- function(x, details = FALSE, ...) {
  print.trunc(x, details, ...)
}

#' @export
print.trunc_exp <- function(x, details = FALSE, ...) {
  print.trunc(x, details, ...)
}

#' @export
print.trunc_gamma <- function(x, details = FALSE, ...) {
  print.trunc(x, details, ...)
}

#' @export
print.trunc_invgamma <- function(x, details = FALSE, ...) {
  print.trunc(x, details, ...)
}

#' @export
print.trunc_invgauss <- function(x, details = FALSE, ...) {
  print.trunc(x, details, ...)
}

#' @export
print.trunc_lognormal <- function(x, details = FALSE, ...) {
  print.trunc(x, details, ...)
}

#' @export
print.trunc_nbinom <- function(x, details = FALSE, ...) {
  print.trunc(x, details, ...)
}

#' @export
print.trunc_normal <- function(x, details = FALSE, ...) {
  print.trunc(x, details, ...)
}

#' @export
print.trunc_poisson <- function(x, details = FALSE, ...) {
  print.trunc(x, details, ...)
}
