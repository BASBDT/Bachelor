.tbc_find_best_num_cut <- function(Xnum, y, control = list()) {
  # --- normalize inputs -------------------------------------------------------
  if (is.data.frame(Xnum)) {
    Xnum <- as.matrix(Xnum)
  }
  if (!is.matrix(Xnum) || !is.numeric(Xnum)) {
    stop("tbc_find_best_num_cut: 'Xnum' must be a numeric matrix/data.frame.", call. = FALSE)
  }
  if (!is.factor(y)) {
    stop("tbc_find_best_num_cut: 'y' must be a factor.", call. = FALSE)
  }

  K <- nlevels(y)
  if (K < 2L) {
    return(NULL)
  }

  # --- resolve class priors ---------------------------------------------------
  # Default: equal priors
  prior <- rep_len(1, K)

  # Prefer high-level name 'prior_probs', fall back to 'prior'
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
    } else if (length(p) == 1L && is.finite(p) && p >= 0) {
      prior <- rep_len(p, K)
    } else {
      warning(
        "tbc_find_best_num_cut: ignoring invalid prior specification; ",
        "using equal class weights."
      )
    }
  }

  # --- optional feature order (1-based indices) -------------------------------
  feature_order <- NULL
  if (!is.null(control$feature_order)) {
    feature_order <- as.integer(control$feature_order)
  }

  # --- C++ call ---------------------------------------------------------------
  res <- tbc_best_numeric_splits(
    X             = Xnum,
    y             = as.integer(y),   # classes 1..K
    priorRatio    = prior,
    feature_order = feature_order    # may be NULL
  )

  # --- guard / mapping --------------------------------------------------------
  if (is.null(res)) return(NULL)
  if (!is.finite(res$bestCrit) ||
      !is.finite(res$bestCut)  ||
      is.null(res$bestVar)     ||
      res$bestVar <= 0L) {
    return(NULL)
  }

  j <- as.integer(res$bestVar)
  nm <- colnames(Xnum)
  feature <- if (!is.null(nm) && length(nm) >= j) nm[[j]] else paste0("X", j)
  naL <- if (!is.null(res$naLeft)) isTRUE(res$naLeft) else TRUE

  list(
    improvement = as.numeric(res$bestCrit),
    cut         = as.numeric(res$bestCut),
    j           = j,
    feature     = feature,
    NA_left     = naL
  )
}
