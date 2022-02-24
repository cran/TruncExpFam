## --##--##--##--##--##--##--##--##--##--##--##--##--
# Evaluating the derivative of log(Gamma(x))
# See https://www.nbi.dk/~polesen/borel/node5.html
# C:  Euler-Mascheroni Constant.
psi <- function(x, k = 10000) {
  # OBS: This function is not used!
  sq <- 1:10000000
  C <- sum(1/sq)-log(max(sq))
  -C + sum(1 / (1:k) - 1 / ((1:k) + x - 1))
}

log.gamma <- function(x, k = 10000) {
  sq <- 1:10000000
  C <- sum(1/sq)-log(max(sq))
  -C * x - log(x) + sum(x / (1:k) - log(1 + x / (1:k)))
}

dpsi.dx <- function(x, k = 10000) {
  # Returns the derivative of the psi function above
  sum((1 / ((0:k) + x))^2)
}
