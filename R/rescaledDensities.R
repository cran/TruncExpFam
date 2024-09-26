rescaledDensities <- function(y, a, b, densFunction, probFunction, ...) {
  dens <- ifelse(
    test = (y < a) | (y > b),
    yes  = 0,
    no   = densFunction(y, ...)
  )
  F.a <- tryCatch(probFunction(a, ...), error = function(e) 0)
  F.b <- ifelse(is.infinite(b), 1, probFunction(b, ...))
  dens <- dens / (F.b - F.a)
  attributes(dens) <- attributes(y)
  dens
}
