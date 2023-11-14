validateFamilyName <- function(family) {
  family <- tolower(family)
  if (!(family %in% valid_distros)) {
    stop(
      "Invalid distribution family. Please choose from the list below:\n",
      paste(valid_distros, collapse = ", ")
    )
  }
}

#' @title Validate family parameters
#' @description Checks if a combination of distribution family and parameters is
#' valid.
#' @param family character with family distribution name
#' @param parms character vector with distribution parameter names
#' @return list telling if family-parm combo is valid + the family name
#' @author Waldir Leoncio
validateFamilyParms <- function(family, parms) {
  matched <- list(family = FALSE, parameters = FALSE)
  families <- grep(family, names(valid_fam_parm))
  for (fam in families) {
    if (any(family == valid_fam_parm[[fam]]$family)) {
      matched$family <- TRUE
      family <- useStandardFamilyName(family)
      parms_text <- paste(parms, collapse = ", ")
      parms_expected <- valid_fam_parm[[fam]]$parms
      if (any(parms == "")) {
        empty_parms <- parms == ""
        parms[empty_parms] <- parms_expected[empty_parms]
        warning(
          "Not all parameters were explicitly specified. ",
          "Please name all function arguments to avoid errors"
        )
      }
      if (all(sort(parms) == sort(parms_expected))) {
        matched$parameters <- TRUE
      }
    }
  }
  if ("parms_expected" %in% ls() && !matched$parameters) {
    parms_expected_text <- paste(unlist(parms_expected), collapse = ", ")
    stop(
      "The {", parms_text, "} ",
      "parameter set does not match the ", family, " family. ",
      "Expected set of parameters: {", parms_expected_text, "}. ",
      "Please change the family to match the expected ",
      "parameters or use a different family."
    )
  }
  return(list(is_valid = all(unlist(matched)), family_name = family))
}

useStandardFamilyName <- function(family) {
  validateFamilyName(family)
  valid_names <- valid_fam_parm[[match(family, valid_distros)]]
  standard_name <- valid_names$family[1]
  return(standard_name)
}
