#' @rdname rtrunc
#' @inheritParams genrtruncClass
#' @export
rtrunc_direct <- function(n, family = "gaussian", parms, a, b, ...) {

  # Validating ---------------------------------------------------------------
  family <-  useStandardFamilyName(tolower(family))
  class(family) <- paste0("trunc_", family)
  validateFamilyName(family)
  validateDomain(family, parms)

  # Determining object class -------------------------------------------------
  class(n) <- genrtruncClass(n, family, names(parms))

  # Dispatching to appropriate sampling funcion ------------------------------
  UseMethod("rtrunc_direct", n)
}

#' @export
rtrunc_direct.normal <- function(n, family, parms, a = -Inf, b = Inf, ...) {
  F_a <- cumDens(a, pnorm, family, parms[["mean"]], parms[["sd"]])
  F_b <- cumDens(b, pnorm, family, parms[["mean"]], parms[["sd"]])
  parms <- c(parms, "n" = n, "a" = a, "b" = b)
  q_T <- truncated_q(
    q_T = qnorm(rescaled_q(n, F_a, F_b), parms[["mean"]], parms[["sd"]]),
    family = family,
    parms = parms
  )
  return(q_T)
}

#' @export
rtrunc_direct.beta <- function(n, family, parms, a = 0, b = 1, ...) {
  F_a <- cumDens(a, pbeta, family, parms[["shape1"]], parms[["shape2"]])
  F_b <- cumDens(b, pbeta, family, parms[["shape1"]], parms[["shape2"]])
  parms <- c(parms, "n" = n, "a" = a, "b" = b)
  q_T <- truncated_q(
    q_T = qbeta(rescaled_q(n, F_a, F_b), parms[["shape1"]], parms[["shape2"]]),
    family = family,
    parms = parms
  )
  return(q_T)
}

#' @export
rtrunc_direct.chisq <- function(n, family, parms, a = 0, b = Inf, ...) {
  F_a <- cumDens(a, pchisq, family, parms[["df"]])
  F_b <- cumDens(b, pchisq, family, parms[["df"]])
  parms <- c(parms, "n" = n, "a" = a, "b" = b)
  q_T <- truncated_q(
    q_T = qchisq(rescaled_q(n, F_a, F_b), parms[["df"]]),
    family = family,
    parms = parms
  )
  return(q_T)
}

#' @export
rtrunc_direct.contbern <- function(n, family, parms, a = 0, b = 1, ...) {
  F_a <- cumDens(a, pcontbern, family, parms[["lambda"]])
  F_b <- cumDens(b, pcontbern, family, parms[["lambda"]])
  parms <- c(parms, "n" = n, "a" = a, "b" = b)
  q_T <- truncated_q(
    q_T = qcontbern(rescaled_q(n, F_a, F_b), parms[["lambda"]]),
    family = family,
    parms = parms
  )
  return(q_T)
}

#' @export
rtrunc_direct.exp <- function(n, family, parms, a = 0, b = Inf, ...) {
  F_a <- cumDens(a, pexp, family, parms[["rate"]])
  F_b <- cumDens(b, pexp, family, parms[["rate"]])
  parms <- c(parms, "n" = n, "a" = a, "b" = b)
  q_T <- truncated_q(
    q_T = qexp(rescaled_q(n, F_a, F_b), parms[["rate"]]),
    family = family,
    parms = parms
  )
  return(q_T)
}

#' @export
rtrunc_direct.gamma <- function(n, family, parms, a = 0, b = Inf, ...) {
  F_a <- cumDens(a, pgamma, family, parms[["shape"]], parms[["rate"]])
  F_b <- cumDens(b, pgamma, family, parms[["shape"]], parms[["rate"]])
  parms <- c(parms, "n" = n, "a" = a, "b" = b)
  if (is.null(parms[["scale"]])) {
    parms[["scale"]] <- 1 / parms[["rate"]]
  } else if (is.null(parms[["rate"]])) {
    parms[["rate"]] <- 1 / parms[["scale"]]
  }
  q_T <- truncated_q(
    q_T = qgamma(rescaled_q(n, F_a, F_b), parms[["shape"]], parms[["rate"]]),
    family = family,
    parms = parms
  )
  return(q_T)
}

#' @export
rtrunc_direct.invgamma <- function(n, family, parms, a = 0, b = Inf, ...) {
  F_a <- cumDens(a, pinvgamma, family, parms[["shape"]], parms[["rate"]])
  F_b <- cumDens(b, pinvgamma, family, parms[["shape"]], parms[["rate"]])
  parms <- c(parms, "n" = n, "a" = a, "b" = b)
  if (is.null(parms[["scale"]])) {
    parms[["scale"]] <- 1 / parms[["rate"]]
  } else if (is.null(parms[["rate"]])) {
    parms[["rate"]] <- 1 / parms[["scale"]]
  }
  q_T <- truncated_q(
    q_T = qinvgamma(rescaled_q(n, F_a, F_b), parms[["shape"]], parms[["rate"]]),
    family = family,
    parms = parms
  )
  return(q_T)
}

#' @export
rtrunc_direct.invgauss <- function(n, family, parms, a = 0, b = Inf, ...) {
  F_a <- cumDens(a, pinvgauss, family, parms[["m"]], parms[["s"]])
  F_b <- cumDens(b, pinvgauss, family, parms[["m"]], parms[["s"]])
  parms <- c(parms, "n" = n, "a" = a, "b" = b)
  q_T <- truncated_q(
    q_T = qinvgauss(rescaled_q(n, F_a, F_b), parms[["m"]], parms[["s"]]),
    family = family,
    parms = parms
  )
  return(q_T)
}

#' @export
rtrunc_direct.lognormal <- function(n, family, parms, a = 0, b = Inf, ...) {
  F_a <- cumDens(a, plnorm, family, parms[["meanlog"]], parms[["sdlog"]])
  F_b <- cumDens(b, plnorm, family, parms[["meanlog"]], parms[["sdlog"]])
  parms <- c(parms, "n" = n, "a" = a, "b" = b)
  q_T <- truncated_q(
    q_T = qlnorm(rescaled_q(n, F_a, F_b), parms[["meanlog"]], parms[["sdlog"]]),
    family = family,
    parms = parms
  )
  return(q_T)
}

#' @export
rtrunc_direct.poisson <- function(n, family, parms, a = 0, b = Inf, ...) {
  F_a <- cumDens(a, ppois, family, parms[["lambda"]])
  F_b <- cumDens(b, ppois, family, parms[["lambda"]])
  # Choose a practical b because a:Inf doesn't work
  practical_b <- ifelse(
    test = b == Inf,
    yes  = qpois(p = 1e-50, lambda = parms[["lambda"]], lower.tail = FALSE),
    no   = b
  )
  parms <- c(parms, "n" = n, "a" = a, "b" = b)
  weights <- dpois(a:practical_b, parms[["lambda"]]) / (F_b - F_a)
  trunc_samp <- sample(a:practical_b, size = n, replace = TRUE, prob = weights)
  f_T <- truncated_q(trunc_samp, family, parms)  # just to add the attributes
  return(f_T)
}

#' @export
rtrunc_direct.binomial <- function(
  n, family, parms, a = 0, b = parms[["size"]], ...
) {
  F_a <- cumDens(a, pbinom, family, parms[["size"]], parms[["prob"]])
  F_b <- cumDens(b, pbinom, family, parms[["size"]], parms[["prob"]])
  parms <- c(parms, "n" = n, "a" = a, "b" = b)
  weights <- dbinom(a:b, parms[["size"]], parms[["prob"]]) / (F_b - F_a)
  trunc_samp <- sample(a:b, size = n, replace = TRUE, prob = weights)
  f_T <- truncated_q(trunc_samp, family, parms)  # just to add the attributes
  return(f_T)
}

#' @export
rtrunc_direct.nbinom <- function(n, family, parms, a = 0, b = Inf, ...) {
  if (is.null(parms[["prob"]])) {
    parms[["prob"]] <- (parms[["size"]]) / (parms[["size"]] + parms[["mu"]])
    parms[["mu"]] <- NULL
  }
  if (is.null(parms[["mu"]])) {
    F_a <- cumDens(a, pnbinom, family, parms[["size"]], parms[["prob"]])
    F_b <- cumDens(b, pnbinom, family, parms[["size"]], parms[["prob"]])
    # Choose a practical b because a:Inf doesn't work
    practical_b <- ifelse(
      test = b == Inf,
      yes  = qnbinom(
        p = 1e-50, parms[["size"]], parms[["prob"]], lower.tail = FALSE
      ),
      no   = b
    )
    d_a_practical_b <- dnbinom(a:practical_b, parms[["size"]], parms[["prob"]])
    weights <- d_a_practical_b / (F_b - F_a)
  }
  trunc_samp <- sample(a:practical_b, size = n, replace = TRUE, prob = weights)
  parms <- c(parms, "n" = n, "a" = a, "b" = b)
  f_T <- truncated_q(trunc_samp, family, parms)  # just to add the attributes
  return(f_T)
}

cumDens <- function(x, probFunction, family, ...) {
  family <- useStandardFamilyName(family)
  support <- valid_fam_parm[[family]][["support"]]
  if (x == min(support)) {
    return(0)
  } else if (x == max(support)) {
    return(1)
  }
  probFunction(x, ...)
}

truncated_q <- function(q_T, family, parms) {
  family <- useStandardFamilyName(family)
  class(q_T) <- paste0("trunc_", family)
  q_T <- attachDistroAttributes(
    sample = q_T,
    family = family,
    parms  = parms
  )
}

rescaled_q <- function(n, F_a, F_b) {
  runif(n) * (F_b - F_a) + F_a
}
