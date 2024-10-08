% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/mlEstimationTruncDist.R
\name{mlEstimationTruncDist}
\alias{mlEstimationTruncDist}
\title{ML Estimation of Distribution Parameters}
\usage{
mlEstimationTruncDist(
  y,
  y.min = attr(y, "truncation_limits")$a,
  y.max = attr(y, "truncation_limits")$b,
  tol = 1e-05,
  max.it = 100,
  delta = 0.33,
  print.iter = 0,
  ny = 100,
  family = NULL,
  ...
)
}
\arguments{
\item{y}{Sequence spanning the domain of the truncated distribution}

\item{y.min}{Lower bound for \code{y}}

\item{y.max}{Upper bound for \code{y}}

\item{tol}{Error tolerance for parameter estimation}

\item{max.it}{Maximum number of iterations}

\item{delta}{Indirectly, the difference between consecutive iterations to
compare with the error tolerance}

\item{print.iter}{Determines the frequency of printing
(i.e., prints every \code{print.iter} iterations)}

\item{ny}{size of intermediate \code{y} range sequence. Higher values yield better
estimations but slower iterations}

\item{family}{distribution family to use}

\item{...}{other parameters passed to subfunctions}
}
\value{
A vector of class \verb{trunc_*} containing the maximum-likelihood
estimation of the underlying distribution * parameters.
}
\description{
ML-estimation of the parameters of the distribution of the
specified family, truncated at y.min and y.max
}
\details{
If \code{print.iter = TRUE}, the function prints the iteration,
the sum of squares of \code{delta.eta.j} (\code{delta.L2}), and the current
parameter estimates. The \code{delta} argument of this function is a factor
in the calculation of \code{delta.eta.j}, which in turn is a factor in the
calculation of \code{delta.L2}.
}
\examples{
sample_size <- 1000
# Normal
sample.norm <- rtrunc(n = sample_size, mean = 2, sd = 1.5, a = -1)
mlEstimationTruncDist(
  sample.norm,
  y.min = -1, max.it = 500, delta = 0.33,
  print.iter = TRUE
)

# Log-Normal
sample.lognorm <- rtrunc(
  n = sample_size, family = "lognormal", meanlog = 2.5, sdlog = 0.5, a = 7
)
ml_lognormal <- mlEstimationTruncDist(
  sample.lognorm,
  y.min = 7, max.it = 500, tol = 1e-10, delta = 0.3,
  print.iter = FALSE
)
ml_lognormal

# Poisson
sample.pois <- rtrunc(
 n = sample_size, lambda = 10, a = 4, family = "Poisson"
)
mlEstimationTruncDist(
  sample.pois,
  y.min = 4, max.it = 500, delta = 0.33,
  print.iter = 5
)

# Gamma
sample.gamma <- rtrunc(
 n = sample_size, shape = 6, rate = 2, a = 2, family = "Gamma"
)
mlEstimationTruncDist(
  sample.gamma,
  y.min = 2, max.it = 1500, delta = 0.3,
  print.iter = 10
)

# Negative binomial
sample.nbinom <- rtruncnbinom(
 sample_size, size = 50, prob = .3, a = 100, b = 120
)
mlEstimationTruncDist(sample.nbinom, r=10)
}
\references{
Inspired by Salvador: Pueyo: "Algorithm for the
maximum likelihood estimation of the parameters of the truncated normal and
lognormal distributions"
}
\author{
René Holst
}
