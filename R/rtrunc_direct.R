#' @rdname rtrunc
#' @param parms list of distribution parameters
#' @export
rtrunc_direct <- function(n, family = "gaussian", parms, ...) {

  # Validating ---------------------------------------------------------------
  family <- tolower(family)
  validateFamilyName(family)

  # Determining object class -------------------------------------------------
  class(n) <- genrtruncClass(n, family, names(parms))

  # Dispatching to appropriate sampling funcion ------------------------------
  UseMethod("rtrunc_direct", n)
}

#' @export
rtrunc_direct.normal <- function(n, family, parms, a = -Inf, b = Inf, ...) {
  F_a <- cumDens(a, pnorm, parms[["mean"]], parms[["sd"]])
  F_b <- cumDens(b, pnorm, parms[["mean"]], parms[["sd"]])
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
  F_a <- cumDens(a, pbeta, parms[["shape1"]], parms[["shape2"]])
  F_b <- cumDens(b, pbeta, parms[["shape1"]], parms[["shape2"]])
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
  F_a <- cumDens(a, pchisq, parms[["df"]])
  F_b <- cumDens(b, pchisq, parms[["df"]])
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
  F_a <- cumDens(a, pcontbern, parms[["lambda"]])
  F_b <- cumDens(b, pcontbern, parms[["lambda"]])
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
  F_a <- cumDens(a, pexp, parms[["rate"]])
  F_b <- cumDens(b, pexp, parms[["rate"]])
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
  F_a <- cumDens(a, pgamma, parms[["shape"]], parms[["rate"]])
  F_b <- cumDens(b, pgamma, parms[["shape"]], parms[["rate"]])
  parms <- c(parms, "n" = n, "a" = a, "b" = b)
  if (is.null(parms[["scale"]])) {
    parms[["scale"]] <- 1 / parms[["rate"]]
  } else if (is.null(parms[["rate"]])) {
    parms[["rate"]] <- 1 / parms[["scale"]]
  } else if (parms[["scale"]] != 1 / parms[["rate"]]) {
    rate <- 1 / parms[["scale"]]
    parms[["rate"]] <- rate
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
  F_a <- cumDens(a, pinvgamma, parms[["shape"]], parms[["rate"]])
  F_b <- cumDens(b, pinvgamma, parms[["shape"]], parms[["rate"]])
  parms <- c(parms, "n" = n, "a" = a, "b" = b)
  if (is.null(parms[["scale"]])) {
    parms[["scale"]] <- 1 / parms[["rate"]]
  } else if (is.null(parms[["rate"]])) {
    parms[["rate"]] <- 1 / parms[["scale"]]
  } else if (parms[["scale"]] != 1 / parms[["rate"]]) {
    rate <- 1 / parms[["scale"]]
    parms[["rate"]] <- rate
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
  F_a <- cumDens(a, pinvgauss, parms[["m"]], parms[["s"]])
  F_b <- cumDens(b, pinvgauss, parms[["m"]], parms[["s"]])
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
  F_a <- cumDens(a, plnorm, parms[["meanlog"]], parms[["sdlog"]])
  F_b <- cumDens(b, plnorm, parms[["meanlog"]], parms[["sdlog"]])
  parms <- c(parms, "n" = n, "a" = a, "b" = b)
  q_T <- truncated_q(
    q_T = qlnorm(rescaled_q(n, F_a, F_b), parms[["meanlog"]], parms[["sdlog"]]),
    family = family,
    parms = parms
  )
  return(q_T)
}

cumDens <- function(x, probFunction, ...) {
  if (x == -Inf || x == 0) {
    return(0)
  } else if (x == Inf || x == 1) {
    return(1)
  } else {
    return(probFunction(x, ...))
  }
}

truncated_q <- function(q_T, family, parms) {
  family <- useStandardFamilyName(family)
  class(q_T) <- paste0("trunc_", family)
  q_T <- attachDistroAttributes(
    sample = q_T,
    family = family,
    parms  = parms
  )
  return(q_T)
}

rescaled_q <- function(n, F_a, F_b) {
  return(runif(n) * (F_b - F_a) + F_a)
}

#' @export
rtrunc_direct.poisson <- function(n, family, parms, a = 0, b = Inf, ...) {
  F_a <- cumDens(a, ppois, parms[["lambda"]])
  F_b <- cumDens(b, ppois, parms[["lambda"]])
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
rtrunc_direct.binomial <- function(n, family, parms, a = 0, b = parms[["size"]], ...) {
  F_a <- cumDens(a, pbinom, parms[["size"]], parms[["prob"]])
  F_b <- cumDens(b, pbinom, parms[["size"]], parms[["prob"]])
  parms <- c(parms, "n" = n, "a" = a, "b" = b)
  weights <- dbinom(a:b, parms[["size"]], parms[["prob"]]) / (F_b - F_a)
  trunc_samp <- sample(a:b, size = n, replace = TRUE, prob = weights)
  f_T <- truncated_q(trunc_samp, family, parms)  # just to add the attributes
  return(f_T)
}

#' @export
rtrunc_direct.nbinom <- function(n, family, parms, a = 0, b = Inf, ...) {
  if (is.null(parms[["mu"]])) {
    F_a <- cumDens(a, pnbinom, parms[["size"]], parms[["prob"]])
    F_b <- cumDens(b, pnbinom, parms[["size"]], parms[["prob"]])
    # Choose a practical b because a:Inf doesn't work
    practical_b <- ifelse(
      test = b == Inf,
      yes  = qnbinom(p = 1e-50, parms[["size"]], parms[["prob"]], lower.tail = FALSE),
      no   = b
    )
    weights <- dnbinom(a:practical_b, parms[["size"]], parms[["prob"]]) / (F_b - F_a)
  } else if (is.null(parms[["prob"]])) {
    F_a <- cumDens(a, pnbinom, parms[["size"]], mu = parms[["mu"]])
    F_b <- cumDens(b, pnbinom, parms[["size"]], mu = parms[["mu"]])
    # Choose a practical b because a:Inf doesn't work
    practical_b <- ifelse(
      test = b == Inf,
      yes  = qnbinom(p = 1e-50, parms[["size"]],  mu = parms[["mu"]], lower.tail = FALSE),
      no   = b
    )
    weights <- dnbinom(a:practical_b, parms[["size"]], mu = parms[["mu"]]) / (F_b - F_a)

  }
  trunc_samp <- sample(a:practical_b, size = n, replace = TRUE, prob = weights)
  parms <- c(parms, "n" = n, "a" = a, "b" = b)
  f_T <- truncated_q(trunc_samp, family, parms)  # just to add the attributes
  return(f_T)
}
