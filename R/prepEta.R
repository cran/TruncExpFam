prepEta <- function(parms, trunc_class) {
  eta <- structure(parms, class = trunc_class)
  if (length(parms) == 1) {
    names(eta) <- "eta"
  } else {
    names(eta) <- paste0("eta", seq_along(parms))
  }
  return(eta)
}
