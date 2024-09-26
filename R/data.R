valid_fam_parm <- list(
  beta = list(
    family = "beta",
    parms = c("shape1", "shape2"),
    support = c(0, 1),
    cont = TRUE
  ),
  binomial = list(
    family = "binomial",
    parms = c("size", "prob"),
    support = c(0, Inf), # Conservative. Upper limit is actually size.
    cont = FALSE
  ),
  chisq = list(
    family = "chisq",
    parms = "df",
    support = c(0, Inf),
    cont = TRUE
  ),
  contbern = list(
    family = "contbern",
    parms = "lambda",
    support = c(0, 1),
    cont = TRUE
  ),
  exp = list(
    family = "exp",
    parms = "rate",
    support = c(0, Inf),
    cont = TRUE
  ),
  gamma = list(
    family = "gamma",
    parms = c("shape", "rate"),
    support = c(0, Inf),
    cont = TRUE
  ),
  gamma = list(
    family = "gamma",
    parms = c("shape", "scale"),
    support = c(0, Inf),
    cont = TRUE
  ),
  invgamma = list(
    family = "invgamma",
    parms = c("shape", "rate"),
    support = c(0, Inf),
    cont = TRUE
  ),
  invgamma = list(
    family = "invgamma",
    parms = c("shape", "scale"),
    support = c(0, Inf),
    cont = TRUE
  ),
  invgauss = list(
    family = "invgauss",
    parms = c("m", "s"),
    support = c(0, Inf),
    cont = TRUE
  ),
  lognormal = list(
    family = "lognormal",
    parms = c("meanlog", "sdlog"),
    support = c(0, Inf),
    cont = TRUE
  ),
  nbinom = list(
    family = "nbinom",
    parms = c("size", "prob"),
    support = c(0, Inf),
    cont = FALSE
  ),
  nbinom = list(
    family = "nbinom",
    parms = c("size", "mu"),
    support = c(0, Inf),
    cont = FALSE
  ),
  normal = list(
    family = c("normal", "gaussian"),
    parms = c("mean", "sd"),
    support = c(-Inf, Inf),
    cont = TRUE
  ),
  gaussian = list(
    family = c("normal", "gaussian"),
    parms = c("mean", "sd"),
    support = c(-Inf, Inf),
    cont = TRUE
  ),
  poisson = list(
    family = "poisson",
    parms = "lambda",
    support = c(0, Inf),
    cont = FALSE
  )
)
valid_distros <- names(valid_fam_parm)
