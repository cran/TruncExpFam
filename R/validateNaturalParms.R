validateNaturalParms <- function(parms) {
  # TODO: develop from scratch or modify validateFamilyParms()?
  if (!is.null(parms) && any(substring(parms, 1, 3) != "eta")) {
    stop("Parameter names must either be empty or start with \"eta\"")
  }
}
