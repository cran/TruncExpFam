attachDistroAttributes <- function(sample, family, parms) {
  if (length(attributes(sample)) == 1) {
    family_attributes <- valid_fam_parm[[family]]
    attr(sample, "parameters") <- parms[family_attributes$parms]
    attr(sample, "truncation_limits") <- parms[c("a", "b")]
    attr(sample, "continuous") <- family_attributes$cont
  }
  return(sample)
}
