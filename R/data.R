valid_fam_parm <- list(
  beta = list(
    family = "beta",
    parms = c("shape1", "shape2"),
    cont = TRUE
  ),
  binomial = list(
    family = "binomial",
    parms = c("size", "prob"),
    cont = FALSE
  ),
  chisq = list(
    family = "chisq",
    parms = "df",
    cont = TRUE
  ),
  contbern = list(
    family = "contbern",
    parms = "lambda",
    cont = TRUE
  ),
  exp = list(
    family = "exp",
    parms = "rate",
    cont = TRUE
  ),
  gamma = list(
    family = "gamma",
    parms = c("shape", "rate"),
    cont = TRUE
  ),
  gamma = list(
    family = "gamma",
    parms = c("shape", "scale"),
    cont = TRUE
  ),
  invgamma = list(
    family = "invgamma",
    parms = c("shape", "rate"),
    cont = TRUE
  ),
  invgamma = list(
    family = "invgamma",
    parms = c("shape", "scale"),
    cont = TRUE
  ),
  invgauss = list(
    family = "invgauss",
    parms = c("m", "s"),
    cont = TRUE
  ),
  lognormal = list(
    family = "lognormal",
    parms = c("meanlog", "sdlog"),
    cont = TRUE
  ),
  nbinom = list(
    family = "nbinom",
    parms = c("size", "prob"),
    cont = FALSE
  ),
  nbinom = list(
    family = "nbinom",
    parms = c("size", "mu"),
    cont = FALSE
  ),
  normal = list(
    family = c("normal", "gaussian"),
    parms = c("mean", "sd"),
    cont = TRUE
  ),
  gaussian = list(
    family = c("normal", "gaussian"),
    parms = c("mean", "sd"),
    cont = TRUE
  ),
  poisson = list(
    family = "poisson",
    parms = "lambda",
    cont = FALSE
  )
)
valid_distros <- names(valid_fam_parm)
