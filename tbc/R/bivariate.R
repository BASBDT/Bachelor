#' Build the internal/canonical label: name_i . a±xx.xx . name_j
#' (kept for backwards-compat / parsing)
#' @keywords internal
.tbc_bivar_label <- function(names, i, j, a, sep = ".") {
  coef_tag <- sprintf("a%0.2f", round(a, 2))
  coef_tag <- sub("^a-", "n", coef_tag)
  paste0(names[i], sep, coef_tag, sep, names[j])
}

#' Build a pretty label for display: "name_i ± xx.xx*name_j"
#' @keywords internal
.tbc_bivar_pretty <- function(names, i, j, a) {
  if (is.na(a)) return(NA_character_)
  s <- if (a >= 0) " + " else " - "
  paste0(names[i], s, sprintf("%0.2f", abs(a)), "*", names[j])
}


#' Find the best bivariate (oblique) split using the C++ backend
#'
#' Thin wrapper around \code{tbc_best_bivariate_split} that:
#' - calls the Rcpp/RcppParallel implementation,
#' - packages indices and the coefficient into \code{DerivedSpec},
#' - returns a structured result consumed by the tree grower,
#' - builds a canonical label using the internal helpers.
#'
#' @param Xnum data.frame or matrix of numeric predictors (no NAs filtered here)
#' @param y factor class vector
#' @param control list produced by \code{tbc_control()}, must contain \code{bivar_coefs}
#' @return list with fields:
#'   \describe{
#'     \item{improvement}{numeric gain/criterion improvement (higher is better)}
#'     \item{cut}{numeric threshold on the projected axis}
#'     \item{i}{primary variable (1-based index)}
#'     \item{j}{secondary variable (1-based index)}
#'     \item{coef}{numeric slope applied to \code{x_j} (alias of \code{a})}
#'     \item{a}{numeric slope used in the linear combination (for internal use)}
#'     \item{is_derived}{logical TRUE}
#'     \item{DerivedSpec}{list(i = ..., j = ..., a = ...)}
#'     \item{DecisionPlaneName}{canonical label built via \code{.tbc_bivar_label}}
#'   }
#' @keywords internal
findCurrentBestBivariateCut <- function(Xnum, y, control) {
  if (!is.factor(y)) {
    stop("findCurrentBestBivariateCut: 'y' must be a factor.", call. = FALSE)
  }

  if (is.data.frame(Xnum)) {
    Xnum <- as.matrix(Xnum)
  }
  if (!is.matrix(Xnum) || !is.numeric(Xnum)) {
    stop("findCurrentBestBivariateCut: 'Xnum' must be a numeric matrix/data.frame.", call. = FALSE)
  }

  if (is.null(control$bivar_coefs) || !length(control$bivar_coefs)) {
    return(NULL)
  }

  coefs <- as.numeric(control$bivar_coefs)
  coefs <- coefs[is.finite(coefs)]
  if (!length(coefs)) return(NULL)

  K <- nlevels(y)

  # default: equal priors
  prior <- rep_len(1, K)

  # allow control to override, preferring names used in the high-level API
  prior_candidate <- NULL
  if (!is.null(control$prior_probs)) {
    prior_candidate <- control$prior_probs
  } else if (!is.null(control$prior)) {
    prior_candidate <- control$prior
  }

  if (!is.null(prior_candidate)) {
    p <- as.numeric(prior_candidate)
    if (length(p) == K && all(is.finite(p)) && all(p >= 0)) {
      prior <- p
    }
  }

  res <- tbc_best_bivariate_split(
    X          = Xnum,
    y          = as.integer(y),
    priorRatio = prior,
    coefs      = coefs
  )

  if (is.null(res)) return(NULL)
  if (!is.finite(res$bestCrit) ||
      is.na(res$bestCut) ||
      is.na(res$i) ||
      is.na(res$j)) {
    return(NULL)
  }

  i    <- as.integer(res$i)
  j    <- as.integer(res$j)
  cut  <- as.numeric(res$bestCut)
  a    <- as.numeric(res$coef)
  naL  <- isTRUE(res$naLeft)

  nm <- colnames(Xnum)
  canon_label <- .tbc_bivar_label(nm, i, j, a)

  list(
    improvement       = as.numeric(res$bestCrit),
    cut               = cut,
    i                 = i,
    j                 = j,
    coef              = a,  # kept for backwards compatibility
    a                 = a,
    NA_left           = naL,
    is_derived        = TRUE,
    DerivedSpec       = list(i = i, j = j, a = a),
    DecisionPlaneName = canon_label
  )
}

#' Internal shim for legacy name used inside tbc_fit
#' @keywords internal
.tbc_find_best_bivar_cut <- function(Xnum, y, control) {
  findCurrentBestBivariateCut(Xnum = Xnum, y = y, control = control)
}
