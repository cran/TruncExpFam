validateSupport <- function(n, ...) {
  UseMethod("validateSupport")
}

#' @method validateSupport trunc_beta
validateSupport.trunc_beta <- function(n, parms, ...) {
  support <- createSupport(n, "[]")
  judgeSupportLimits(n, parms, support)
}

#' @method validateSupport trunc_binomial
validateSupport.trunc_binomial <- function(n, parms, nsize = parms$size, ...) {
  support <- createSupport(n, "{}")
  support[["u"]] <- nsize
  judgeSupportLimits(n, parms, support, FALSE)
}

#' @method validateSupport trunc_chisq
validateSupport.trunc_chisq <- function(n, parms, ...) {
  if (is.null(parms) || parms$df > 1) {
    support <- createSupport(n, "[)")
  } else {
    support <- createSupport(n, "()")
  }
  judgeSupportLimits(n, parms, support)
}

#' @method validateSupport trunc_contbern
validateSupport.trunc_contbern <- function(n, parms, ...) {
  support <- createSupport(n, "[]")
  judgeSupportLimits(n, parms, support)
}

#' @method validateSupport trunc_exp
validateSupport.trunc_exp <- function(n, parms, ...) {
  support <- createSupport(n, "[)")
  judgeSupportLimits(n, parms, support)
}

#' @method validateSupport trunc_gamma
validateSupport.trunc_gamma <- function(n, parms, ...) {
  support <- createSupport(n, "()")
  judgeSupportLimits(n, parms, support)
}

#' @method validateSupport trunc_invgamma
validateSupport.trunc_invgamma <- function(n, parms, ...) {
  support <- createSupport(n, "()")
  judgeSupportLimits(n, parms, support)
}

#' @method validateSupport trunc_invgauss
validateSupport.trunc_invgauss <- function(n, parms, ...) {
  support <- createSupport(n, "()")
  judgeSupportLimits(n, parms, support)
}

#' @method validateSupport trunc_lognormal
validateSupport.trunc_lognormal <- function(n, parms, ...) {
  support <- createSupport(n, "()")
  judgeSupportLimits(n, parms, support)
}

#' @method validateSupport trunc_nbinom
validateSupport.trunc_nbinom <- function(n, parms, ...) {
  support <- createSupport(n, "{}")
  judgeSupportLimits(n, parms, support, FALSE)
}

#' @method validateSupport trunc_normal
validateSupport.trunc_normal <- function(n, parms, ...) {
  support <- createSupport(n, "()")
  judgeSupportLimits(n, parms, support, no_complex = TRUE)
}

#' @method validateSupport trunc_poisson
validateSupport.trunc_poisson <- function(n, parms, ...) {
  support <- createSupport(n, "{}")
  judgeSupportLimits(n, parms, support, FALSE)
}

createSupport <- function(n, inclusion_brackets) {
  # This function outputs a list of numbers and texts related to the support
  # of a truncated distribution. It does not evaluate or validate anything,
  # it just blindly builds the output.
  family <- gsub("trunc_", "", class(n))
  support <- valid_fam_parm[[family]][["support"]]
  out <- list(l = min(support), u = max(support), txt = vector("character"))
  split_brackets <- strsplit(inclusion_brackets, "")
  for (i in seq_along(split_brackets)) {
    lower_symbol <- split_brackets[[i]][1]
    upper_symbol <- split_brackets[[i]][2]
      if (lower_symbol == "(") {
        new_txt <- paste0(lower_symbol, out$l, ", ", out$u, upper_symbol)
      } else {
        if (is.infinite(out$u)) {
          new_txt <- paste0(lower_symbol, out$l, ", ...}")
        } else {
          new_txt <- paste0(lower_symbol, out$l, ", ..., ", out$u, upper_symbol)
        }
      }
    out$txt <- append(out$txt, new_txt)
  }
  out$txt <- paste(out$txt, collapse = " or ")
  out
}

judgeSupportLimits <- function(
  data, parms, support, cont = TRUE, no_complex = FALSE
) {
  # Data circuit breaker =======================================================
  operator_1 <- strsplit(support[["txt"]], "")[[1]][1]
  operator_2 <- strsplit(support[["txt"]], "")[[1]][nchar(support[["txt"]])]
  if (operator_1 %in% c("{", "[")) {
    operator_1 <- get("<")
  } else {
    operator_1 <- get("<=")
  }
  if (operator_2 %in% c("}", "]")) {
    operator_2 <- get(">")
  } else {
    operator_2 <- get(">=")
  }
  if (any(operator_1(data, support$l)) || any(operator_2(data, support$u))) {
    stop("Sample contains values outside of support ", support$txt)
  }

  # Truncation limits circuit breaker ==========================================
  if (!is.null(parms$a) && !is.null(parms$b)) {
    # Treating complex numbers =================================================
    if (no_complex && (is.complex(parms$a) || is.complex(parms$b))) {
      stop("Truncation limits may not contain complex numbers")
    }

    # Treating edge cases ======================================================
    split_brackets <- strsplit(support$txt, "")[[1]]
    if (!is.null(parms$a) && !is.null(parms$b)) {
      if (cont) {
        cond_al <- parms$a < support$l
        cond_au <- parms$a >= support$u # would leave no margin for sampling
        cond_bl <- parms$b <= support$l # would leave no margin for sampling
        cond_bu <- parms$b > support$u
      } else {
        cond_al <- parms$a < support$l
        cond_au <- parms$a > support$u
        cond_bl <- parms$b < support$l
        cond_bu <- parms$b > support$u
      }
    }

    # Judging support limits ===================================================
    if (parms$a == parms$b) {
      stop("Identical truncation limits: a = b = ", parms$a)
    } else if (cond_au || cond_bl) {
      stop(
        "Truncation limits {", parms$a, ", ", parms$b, "} must be a subset of ",
        support$txt
      )
    } else if (cond_al || cond_bu) {
      warning(
        "Truncation limits {", parms$a, ", ", parms$b, "} are not a subset of ",
        support$txt
      )
    } else if (parms$b <= parms$a) {
      stop("Upper truncation limit (b) must be higher than lower limit (a)")
    }
  }

}
