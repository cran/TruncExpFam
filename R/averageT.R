#' @title Averages out the sufficient statistics T(y)
#' @description Takes a vector of values and returns the column average of their
#' sufficient statistic (determined by their class)
#' @param y vector of values
#' @return A vector with the average of the sufficient statistics
averageT <- function(y) {
  Ty <- sufficientT(y)
  if (is.null(dim(Ty))) {
    out <- mean(Ty)
  } else {
    out <- colMeans(Ty)
  }
  return(out)
}
