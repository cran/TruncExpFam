#' @title Cumulative Distribution Function
#' @description Calculates the cumulative probability for a given truncated
#' distribution
#' @param q vector of quantiles
#' @param family distribution family to use
#' @param lower.tail logical; if `TRUE`, probabilities are
#' \eqn{P(X <= x)}{P(X \leq x)} otherwise, \eqn{P(X > x)}
#' @param log.p logical; if `TRUE`, probabilities p are given as `log(p)`
#' @param ... \emph{named} distribution parameters and/or truncation limits
#' (`a`, `b`)
#' @export
#' @return The cumulative probability of y.
#' @examples
#' ptrunc(0)
#' ptrunc(6, family = "gaussian", mean = 5, sd = 10, b = 7)
#' pnorm(6, mean = 5, sd = 10) # for comparison
ptrunc <- function(q, family, ..., lower.tail = TRUE, log.p = FALSE) {
  # This is a pseudo-generic function to pre-process arguments and call the
  # actual generic, ptrunc.generic().

  # Validating ---------------------------------------------------------------
  if (missing(family)) family <- "gaussian"
  family <- tolower(family)
  validateFamilyName(family)

  # Reclassifying q and dispatching -------------------------------------------
  class(q) <- genrtruncClass(q, family, NULL)
  prob <- ptrunc.generic(q, ..., lower.tail = lower.tail, log.p = log.p)
  return(unclass(prob))
}

ptrunc.generic <- function(q, ..., lower.tail, log.p) {
  UseMethod("ptrunc", q)
}

ptrunc.normal <- function(
  q, mean = 0, sd = 1, a = -Inf, b = Inf, ..., lower.tail = TRUE, log.p = FALSE
  ) {
  validate_q_a_b(q, a, b)
  p_q <- pnorm(q, mean, sd, lower.tail = TRUE, log.p)
  p_a <- pnorm(a, mean, sd, lower.tail = TRUE, log.p)
  p_b <- pnorm(b, mean, sd, lower.tail = TRUE, log.p)
  return(truncated_p(p_q, p_a, p_b, lower.tail, log.p))
}

#' @export
#' @rdname ptrunc
#' @inheritParams rtrunc
ptruncnorm <- ptrunc.normal

ptrunc.beta <- function(
  q, shape1, shape2, a = 0, b = 1, ..., lower.tail = TRUE, log.p = FALSE
) {
  validate_q_a_b(q, a, b)
  p_q <- pbeta(q, shape1, shape2, ncp = 0, lower.tail = TRUE, log.p)
  p_a <- pbeta(a, shape1, shape2, ncp = 0, lower.tail = TRUE, log.p)
  p_b <- pbeta(b, shape1, shape2, ncp = 0, lower.tail = TRUE, log.p)
  return(truncated_p(p_q, p_a, p_b, lower.tail, log.p))
}

#' @export
#' @rdname ptrunc
#' @inheritParams rtrunc
ptruncbeta <- ptrunc.beta

ptrunc.binomial <- function(
  q, size, prob, a = 0, b = size, ..., lower.tail = TRUE, log.p = FALSE
) {
  validate_q_a_b(q, a, b)
  p_q <- pbinom(q, size, prob, lower.tail = TRUE, log.p)
  p_a <- pbinom(a - 1L, size, prob, lower.tail = TRUE, log.p)
  p_b <- pbinom(b, size, prob, lower.tail = TRUE, log.p)
  return(truncated_p(p_q, p_a, p_b, lower.tail, log.p))
}

#' @export
#' @rdname ptrunc
#' @inheritParams rtrunc
ptruncbinom <- ptrunc.binomial

ptrunc.poisson <- function(
  q, lambda, a = 0, b = Inf, ..., lower.tail = TRUE, log.p = FALSE
) {
  validate_q_a_b(q, a, b)
  p_q <- ppois(q, lambda, lower.tail = TRUE, log.p)
  p_a <- ppois(a - 1L, lambda, lower.tail = TRUE, log.p)
  p_b <- ppois(b, lambda, lower.tail = TRUE, log.p)
  return(truncated_p(p_q, p_a, p_b, lower.tail, log.p))
}

#' @export
#' @rdname ptrunc
#' @inheritParams rtrunc
ptruncpois <- ptrunc.poisson

ptrunc.chisq <- function(
  q, df, a = 0, b = Inf, ..., lower.tail = TRUE, log.p = FALSE
) {
  validate_q_a_b(q, a, b)
  p_q <- pchisq(q, df, ncp = 0, lower.tail = TRUE, log.p)
  p_a <- pchisq(a - 1L, df, ncp = 0, lower.tail = TRUE, log.p)
  p_b <- pchisq(b, df, ncp = 0, lower.tail = TRUE, log.p)
  return(truncated_p(p_q, p_a, p_b, lower.tail, log.p))
}

#' @export
#' @rdname ptrunc
#' @inheritParams rtrunc
ptruncchisq <- ptrunc.chisq

ptrunc.contbern <- function(q, lambda, a = 0, b = 1, ...) {
  validate_q_a_b(q, a, b)
  p_q <- pcontbern(q, lambda)
  p_a <- pcontbern(a, lambda)
  p_b <- pcontbern(b, lambda)
  return(truncated_p(p_q, p_a, p_b, lower.tail = TRUE, log.p = FALSE))
}

#' @export
#' @rdname ptrunc
#' @inheritParams rtrunc
ptrunccontbern <- ptrunc.contbern

ptrunc.exp <- function(
  q, rate = 1, a = 0, b = Inf, ..., lower.tail = TRUE, log.p = FALSE
) {
  validate_q_a_b(q, a, b)
  p_q <- pexp(q, rate, lower.tail = TRUE, log.p)
  p_a <- pexp(a, rate, lower.tail = TRUE, log.p)
  p_b <- pexp(b, rate, lower.tail = TRUE, log.p)
  return(truncated_p(p_q, p_a, p_b, lower.tail, log.p))
}

#' @export
#' @rdname ptrunc
#' @inheritParams rtrunc
ptruncexp <- ptrunc.exp

ptrunc.gamma <- function(
  q, shape, rate = 1, scale = 1 / rate, a = 0, b = Inf, ..., lower.tail = TRUE,
  log.p = FALSE
) {
  validate_q_a_b(q, a, b)
  if (!missing(rate) && !missing(scale)) {
    stop("specify 'rate' or 'scale' but not both")
  }
  p_q <- pgamma(q, shape, scale = scale, lower.tail = TRUE, log.p = log.p)
  p_a <- pgamma(a, shape, scale = scale, lower.tail = TRUE, log.p = log.p)
  p_b <- pgamma(b, shape, scale = scale, lower.tail = TRUE, log.p = log.p)
  return(truncated_p(p_q, p_a, p_b, lower.tail, log.p))
}

#' @export
#' @rdname ptrunc
#' @inheritParams rtrunc
ptruncgamma <- ptrunc.gamma

ptrunc.invgamma <- function(
  q, shape, rate = 1, scale = 1 / rate, a = 0, b = Inf, ..., lower.tail = TRUE,
  log.p = FALSE
) {
  validate_q_a_b(q, a, b)
  if (!missing(rate) && !missing(scale)) {
    stop("specify 'rate' or 'scale' but not both")
  }
  p_q <- pinvgamma(q, shape, scale = scale, lower.tail = TRUE, log.p = log.p)
  p_a <- pinvgamma(a, shape, scale = scale, lower.tail = TRUE, log.p = log.p)
  p_b <- pinvgamma(b, shape, scale = scale, lower.tail = TRUE, log.p = log.p)
  return(truncated_p(p_q, p_a, p_b, lower.tail, log.p))
}

#' @export
#' @rdname ptrunc
#' @inheritParams rtrunc
ptruncinvgamma <- ptrunc.invgamma

ptrunc.invgauss <- function(q, m, s, a = 0, b = Inf, ...) {
  validate_q_a_b(q, a, b)
  p_q <- pinvgauss(q, m, s)
  p_a <- ifelse(a == 0, 0, pinvgauss(a, m, s))
  p_b <- ifelse(b == Inf, 1, pinvgauss(b, m, s))
  return(truncated_p(p_q, p_a, p_b, lower.tail = TRUE, log.p = FALSE))
}

#' @export
#' @rdname ptrunc
#' @inheritParams rtrunc
ptruncinvgauss <- ptrunc.invgauss

ptrunc.lognormal <- function(
  q, meanlog = 0, sdlog = 1, a = 0, b = Inf, ..., lower.tail = TRUE,
  log.p = FALSE
) {
  validate_q_a_b(q, a, b)
  p_q <- plnorm(q, meanlog, sdlog, lower.tail = TRUE, log.p)
  p_a <- plnorm(a, meanlog, sdlog, lower.tail = TRUE, log.p)
  p_b <- plnorm(b, meanlog, sdlog, lower.tail = TRUE, log.p)
  return(truncated_p(p_q, p_a, p_b, lower.tail, log.p))
}

#' @export
#' @rdname ptrunc
#' @inheritParams rtrunc
ptrunclnorm <- ptrunc.lognormal

ptrunc.nbinom <- function(
  q, size, prob, mu, a = 0, b = Inf, ..., lower.tail = TRUE, log.p = FALSE
) {
  if (missing(prob)) {
    prob <- size / (size + mu) # from help("pnbinom")
    mu <- NULL
  }
  validate_q_a_b(q, a, b)
  p_q <- pnbinom(q, size, prob, lower.tail = TRUE, log.p = log.p)
  p_a <- pnbinom(a - 1L, size, prob, lower.tail = TRUE, log.p = log.p)
  p_b <- pnbinom(b, size, prob, lower.tail = TRUE, log.p = log.p)
  return(truncated_p(p_q, p_a, p_b, lower.tail, log.p))
}

#' @export
#' @rdname ptrunc
#' @inheritParams rtrunc
ptruncnbinom <- ptrunc.nbinom

truncated_p <- function(p_q, p_a, p_b, lower.tail = TRUE, log.p = FALSE) {
  # Usual cases --------------------------------------------------------------
  if (log.p) {
    p <- log((exp(p_q) - exp(p_a)) / (exp(p_b) - exp(p_a)))
    if (!lower.tail) {
      p <- log(1 - exp(p))
    }
  } else {
    p <- (p_q - p_a) / (p_b - p_a)
    if (!lower.tail) {
      p <- 1 - p
    }
  }
  return(p)
}

validate_q_a_b <- function(q, a, b) {
  if (a > b) stop("a must be <= b")
  if (any(q < a) || any(q > b)) stop("q must be in [a, b]")
}
