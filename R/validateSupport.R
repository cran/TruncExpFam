validateSupport <- function(n, ...) {
  UseMethod("validateSupport")
}

validateSupport.trunc_beta <- function(n, parms, ...) {
  support <- createSupport(0, 1, "()")
  judgeSupportLimits(parms, support)
}

validateSupport.trunc_binomial <- function(n, parms, ...) {
  support <- createSupport(0, parms$size, "{}")
  judgeSupportLimits(parms, support, FALSE)
}

validateSupport.trunc_chisq <- function(n, parms, ...) {
  support <- createSupport(0, Inf, "())") # Wikipedia uses [0, Inf) for df > 1
  judgeSupportLimits(parms, support)
}

validateSupport.trunc_contbern <- function(n, parms, ...) {
  support <- createSupport(0, 1, "[]")
  judgeSupportLimits(parms, support)
}

validateSupport.trunc_exp <- function(n, parms, ...) {
  support <- createSupport(0, Inf, "()") # Wikipedia says [0, Inf]; R uses (0,Inf)
  judgeSupportLimits(parms, support)
}

validateSupport.trunc_gamma <- function(n, parms, ...) {
  support <- createSupport(0, Inf, "()")
  judgeSupportLimits(parms, support)
}

validateSupport.trunc_invgamma <- function(n, parms, ...) {
  support <- createSupport(0, Inf, "()")
  judgeSupportLimits(parms, support)
}

validateSupport.trunc_invgauss <- function(n, parms, ...) {
  support <- createSupport(0, Inf, "()")
  judgeSupportLimits(parms, support)
}

validateSupport.trunc_lognormal <- function(n, parms, ...) {
  support <- createSupport(0, Inf, "()")
  judgeSupportLimits(parms, support)
}

validateSupport.trunc_nbinom <- function(n, parms, ...) {
  support <- createSupport(0, Inf, "{}")
  judgeSupportLimits(parms, support, FALSE)
}

validateSupport.trunc_normal <- function(n, parms, ...) {
  support <- createSupport(-Inf, Inf, "()")
  judgeSupportLimits(parms, support, no_complex = TRUE)
}

validateSupport.trunc_poisson <- function(n, parms, ...) {
  support <- createSupport(0, Inf, "{}")
  judgeSupportLimits(parms, support, FALSE)
}

createSupport <- function(lower, upper, inclusion_brackets) {
  out <- list(
    l = lower,
    u = upper,
    txt = vector("character")
  )
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
  return(out)
}

judgeSupportLimits <- function(parms, support, cont = TRUE, no_complex = FALSE) {
  # Complex numbers circuit breaker ============================================
  if (no_complex & (is.complex(parms$a) | is.complex(parms$b))) {
    stop("Truncation limits may not contain complex numbers")
  }

  # Treating edge cases ========================================================
  split_brackets <- strsplit(support$txt, "")[[1]]
  include_l <- split_brackets[1] %in% c("[", "{")
  include_u <- split_brackets[length(split_brackets)] %in% c("]", "}")
  if (cont) {
    cond_au <- parms$a >= support$u
    cond_bl <- parms$b <= support$l
    cond_al <- parms$a < support$l
    cond_bu <- parms$b > support$u
  } else {
    if (include_l) {
      cond_al <- parms$a < support$l
      cond_bl <- parms$b < support$l
    } else {
      cond_al <- parms$a <= support$l
      cond_bl <- parms$b <= support$l
    }
    if (include_u) {
      cond_au <- parms$a > support$u
      cond_bu <- parms$b > support$u
    } else {
      cond_au <- parms$a >= support$u
      cond_bu <- (parms$b >= support$u) & (parms$b != Inf)
    }
  }

  # Judging suppor limits ======================================================
  if (parms$a == parms$b) {
    stop("Identical truncation limits: a = b = ", parms$a)
  } else if (cond_au | cond_bl) {
    stop(
      "Truncation limits {", parms$a, ", ", parms$b, "} must be a subset of ",
      support$txt
    )
  } else if (cond_al | cond_bu) {
    warning(
      "Truncation limits {", parms$a, ", ", parms$b, "} are not a subset of ",
      support$txt
    )
  } else if (parms$b <= parms$a) {
    stop("Upper truncation limit (b) must be higher than lower limit (a)")
  }
}
