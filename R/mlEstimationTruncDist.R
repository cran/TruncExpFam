#' @title ML Estimation of Distribution Parameters
#' @description ML-estimation of the parameters of the distribution of the specified family, truncated at y.min and y.max
#' @param y Sequence spanning the domain of the truncated distribution
#' @param y.min Lower bound for y
#' @param y.max Upper bound for y
#' @param tol Error tolerance for parameter estimation
#' @param delta Indirectly, the difference between consecutive iterations to compare with the error tolerance
#' @param max.it Maximum number of iterations
#' @param print.iter Determines the frequency of printing (i.e., prints every \code{print.iter} iterations)
#' @param ny size of intermediate y range sequence. Higher values yield better estimations but slower iterations
#' @param ... other parameters passed to subfunctions
#' @note `print.iter` can be `TRUE`, `FALSE` or an integer indicating an interval for printing every `X` iterations.
#' @references Inspired by Salvador: Pueyo: "Algorithm for the maximum likelihood estimation of the parameters of the truncated normal and lognormal distributions"
#' @author René Holst
#' @importFrom stats dbinom dgamma dlnorm dnorm dpois pbinom pgamma plnorm pnorm ppois rbinom rgamma rlnorm rnorm rpois var
#' @importFrom methods is
#' @examples
#' sample_size <- 1000
#' # Normal
#' sample.norm <- rtrunc(n = sample_size, mean = 2, sd = 1.5, a = -1)
#' mlEstimationTruncDist(
#'   sample.norm,
#'   y.min = -1, max.it = 500, delta = 0.33,
#'   print.iter = TRUE
#' )
#'
#' # Log-Normal
#' sample.lognorm <- rtrunc(
#'   n = sample_size, family = "lognormal", meanlog = 2.5, sdlog = 0.5, a = 7
#' )
#' ml_lognormal <- mlEstimationTruncDist(
#'   sample.lognorm,
#'   y.min = 7, max.it = 500, tol = 1e-10, delta = 0.3,
#'   print.iter = FALSE
#' )
#' ml_lognormal
#'
#' # Poisson
#' sample.pois <- rtrunc(n = sample_size, lambda = 10, a = 4, family = "Poisson")
#' mlEstimationTruncDist(
#'   sample.pois,
#'   y.min = 4, max.it = 500, delta = 0.33,
#'   print.iter = 5
#' )
#'
#' # Gamma
#' sample.gamma <- rtrunc(n = sample_size, shape = 6, rate = 2, a = 2, family = "Gamma")
#' mlEstimationTruncDist(
#'   sample.gamma,
#'   y.min = 2, max.it = 1500, delta = 0.3,
#'   print.iter = 10
#' )
#'
#' # Negative binomial
#' sample.nbinom <- rtruncnbinom(sample_size, size = 50, prob = .3, a = 100, b = 120)
#' mlEstimationTruncDist(sample.nbinom, r=10)
#' @export
#' @return A vector of class \code{trunc_*} containing the maximum-likelihood
#' estmation of the underlying distribution * parameters.
mlEstimationTruncDist <- function(y, y.min = attr(y, "truncation_limits")$a,
  y.max = attr(y, "truncation_limits")$b, tol = 1e-5, max.it = 25,
  delta = 0.33, print.iter = 0, ny = 100, ...
) {
  # Some initialisations
  if (as.numeric(print.iter) > 0) {
    distro_name <- gsub("trunc_", "", class(y))
    message("Estimating parameters for the ", distro_name, " distribution")
  }
  T.avg <- averageT(y)
  eta.j <- parameters2natural(init.parms(y))
  y.seq <- getYseq(y, y.min, y.max, ny) # y-values to calculate expectations
  it <- 0
  delta.L2 <- 10000 # sum of squares of individual delta.eta.j (see below)
  # Now iterate
  while ((delta.L2 > tol) & (it < max.it)) {
    parm.j <- natural2parameters(eta.j)
    T.minus.E.T <- getTminusET(
      eta.j, y.seq, y.min, y.max, attr(y, "continuous"), T.avg
    )
    grad.E.T.inv <- getGradETinv(eta.j, ...) # p x p
    delta.eta.j.plus.1 <- delta * grad.E.T.inv %*% T.minus.E.T
    eta.j <- eta.j + delta.eta.j.plus.1
    delta.L2 <- sum(delta.eta.j.plus.1^2)
    it <- it + 1
    if (print.iter) {
      if (it %% as.numeric(print.iter) == 0) {
        cat(
          "it: ", it, "delta: ", delta.L2, " - parm: ",
          round(parm.j, 3), "\n"
        )
      }
    }
  }
  parm <- natural2parameters(eta.j)
  return(parm)
}

getTminusET <- function(eta, y.seq, y.min, y.max, cont.dist, T.avg) {
  # Calculates T.bar-E(T|eta_j) by numerical integration
  delta.y <- y.seq[2] - y.seq[1] # step length, length(y.seq)=L
  trunc.density <- dtrunc(y.seq, eta, y.min, y.max) # L vector
  T.f <- sufficientT(y.seq) * trunc.density # L x p matrix
  if (length(eta) > 1) {
    E.T.j <- delta.y * apply(T.f, 2, sum) # 1 x p
    if (cont.dist) {
      E.T.j <- E.T.j - delta.y * 0.5 * (T.f[1, ] + T.f[length(y.seq), ])
    }
  } else {
    E.T.j <- delta.y * sum(T.f)
    if (cont.dist) {
      E.T.j <- E.T.j - delta.y * 0.5 * (T.f[1] + T.f[length(y.seq)])
    }
  }
  T.bar.minus.E.T.j <- T.avg - E.T.j # 1 x p
  return(T.bar.minus.E.T.j)
}
