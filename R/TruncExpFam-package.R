#' @title Truncated Exponential Family
#' @description TruncExpFam is an R package to handle truncated members from the
#' exponential family.
#' @docType package
#' @name TruncExpFam
#' @details This package offers truncated counterparts of the density-,
#' distribution-, quantile- and sampling-functions for a broad range of
#' distributions from the exponential family, as implemented in the [stats]
#' package.
#'
#' The package also provides functions for estimating the parameters of the
#' distributions from data, given the truncation limits.
#'
#' For more info, please check [rtrunc()], [dtrunc()] and [print.trunc()].
#' Counterparts for density and probability functions are on the roadmap for
#' a future release.
#'
#' @section Supported distributions:
#' \itemize{
#'  \item Beta
#'  \item Binomial
#'  \item Chi-Square
#'  \item Continuous Bernoulli
#'  \item Exponential
#'  \item Gamma
#'  \item Inverse Gamma
#'  \item Inverse Gaussian
#'  \item Log-normal
#'  \item Negative Binomial
#'  \item Normal
#'  \item Poisson
#' }
#' @note Found a bug? Want to suggest a feature? Contribute to the scientific
#' and open source communities by opening an issue on our home page.
#' Check the "BugReports" field on `packageDescription("TruncExpFam")` for
#' the URL.
#' @importFrom methods new is
#' @importFrom stats dbinom dgamma dlnorm dnorm dpois pbinom pgamma plnorm pnorm
#' ppois rbinom rgamma rlnorm rnorm rpois var dbeta pbeta rbeta dchisq pchisq
#' qexp qgamma rchisq runif rexp dexp pexp dnbinom pnbinom rnbinom qnorm qbeta
#' qchisq qlnorm qpois qbinom qnbinom
#' @importFrom invgamma rinvgamma dinvgamma pinvgamma qinvgamma
#' @importFrom rmutil rinvgauss dinvgauss pinvgauss qinvgauss
"_PACKAGE"
