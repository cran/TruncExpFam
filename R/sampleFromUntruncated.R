# this is an internal function used by each rtrunc* function to generate initial
# values
sampleFromTruncated <- function(parms) {

  # Generating empty structure with correct class ------------------------------
  className <- class(parms$n)
  family <- gsub("trunc_", "", className)
  y <- structure(numeric(0), class = className)

  # Checking domain and parmeters ----------------------------------------------
  validateDomain(y, parms)
  validateSupport(y, parms)
  if (family %in% c("gamma", "invgamma")) {
    parms$rate <- NULL
  } else if (family == "nbinom") {
    if (parms$mu == "") {
      parms$mu <- NULL
    } else if (parms$prob == "") {
      parms$prob <- NULL
      class(parms$n) <- "trunc_nbinom_mu"
    }
  }
  common_parms <- c("a", "b", "n")
  validateFamilyParms(family, names(parms)[!(names(parms) %in% common_parms)])

  # Sampling from untruncated distro -------------------------------------------
  y <- with(parms, switch(class(n),
    "trunc_beta" = rbeta(n, shape1, shape2),
    "trunc_binomial" = rbinom(n, size, prob),
    "trunc_chisq" = rchisq(n, df),
    "trunc_contbern" = rcontbern(n, lambda),
    "trunc_exp" = rexp(n, rate),
    "trunc_gamma" = rgamma(n, shape = shape, scale = scale),
    "trunc_invgamma" = rinvgamma(n, shape = shape, scale = scale),
    "trunc_invgauss" = rinvgauss(n, m, s),
    "trunc_lognormal" = rlnorm(n, meanlog, sdlog),
    "trunc_nbinom" = rnbinom(n, size, prob),
    "trunc_nbinom_mu" = rnbinom(n, size, mu = mu),
    "trunc_normal" = rnorm(n, mean, sd),
    "trunc_poisson" = rpois(n, lambda)
  ))
  y <- y[y >= parms$a]
  y <- y[y <= parms$b]

  extra_n <- 1 # to generate extra observations to complete n from input
  while (length(y) != parms$n) {
    new_y <- with(parms, switch(class(n),
      "trunc_beta" = rbeta(extra_n, shape1, shape2),
      "trunc_binomial" = rbinom(extra_n, size, prob),
      "trunc_chisq" = rchisq(extra_n, df),
      "trunc_contbern" = rcontbern(extra_n, lambda),
      "trunc_exp" = rexp(extra_n, rate),
      "trunc_gamma" = rgamma(extra_n, shape = shape, scale = scale),
      "trunc_invgamma" = rinvgamma(extra_n, shape = shape, scale = scale),
      "trunc_invgauss" = rinvgauss(extra_n, m, s),
      "trunc_lognormal" = rlnorm(extra_n, meanlog, sdlog),
      "trunc_nbinom" = rnbinom(extra_n, size, prob),
      "trunc_nbinom_mu" = rnbinom(extra_n, size, mu = mu),
      "trunc_normal" = rnorm(extra_n, mean, sd),
      "trunc_poisson" = rpois(extra_n, lambda)
    ))

    # Truncating ---------------------------------------------------------------
    new_y <- new_y[new_y >= parms$a]
    new_y <- new_y[new_y <= parms$b]

    # Assembling ---------------------------------------------------------------
    y <- append(y, new_y)
  }

  # Attaching attributes -------------------------------------------------------
  class(y) <- className
  y <- attachDistroAttributes(y, family, parms)
  return(y)
}
