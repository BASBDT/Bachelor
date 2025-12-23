#' Split Predictors into Numeric and Categorical Matrices (Internal)
#'
#' Internally splits a dataset into numeric and categorical matrices based on
#' type detection and optional overrides.
#'
#' @param data A matrix or data.frame.
#' @param feature_not_numeric Optional mask (logical/integer/character) indicating categorical columns.
#'
#' @return A list with \code{numeric_data} and \code{categoric_data} matrices.
#' @keywords internal
.tbc_split_numeric_categoric <- function(data, feature_not_numeric = NULL) {

  # Normalize Input
  if (is.data.frame(data)) {
    df <- as.data.frame(data, stringsAsFactors = FALSE, optional = TRUE)
  } else if (is.matrix(data)) {
    df <- as.data.frame(data, stringsAsFactors = FALSE)
  } else {
    stop(".tbc_split_numeric_categoric: 'data' must be a matrix or data.frame.", call. = FALSE)
  }

  n_rows <- nrow(df)
  n_cols <- ncol(df)

  # Ensure column names exist
  if (n_cols > 0L && is.null(colnames(df))) {
    colnames(df) <- paste0("V", seq_len(n_cols))
  }

  # Helper to resolve categorical mask
  resolve_mask <- function(dat, idx) {
    d <- ncol(dat)
    nm <- colnames(dat)
    if (is.null(idx)) return(NULL)

    if (is.logical(idx)) {
      if (length(idx) != d) stop("feature_not_numeric (logical) must have length = ncol(data).", call. = FALSE)
      return(idx)
    }

    if (is.numeric(idx) || is.integer(idx)) {
      mask <- rep(FALSE, d)
      ii <- as.integer(idx)
      if (any(!is.finite(ii))) stop("feature_not_numeric contains non-finite indices.", call. = FALSE)
      if (any(ii < 1L | ii > d)) stop("feature_not_numeric indices out of range.", call. = FALSE)
      mask[ii] <- TRUE
      return(mask)
    }

    if (is.character(idx)) {
      miss <- setdiff(idx, nm)
      if (length(miss)) stop("feature_not_numeric names not found in data.", call. = FALSE)
      return(nm %in% idx)
    }

    stop("Unsupported type for feature_not_numeric.", call. = FALSE)
  }

  cat_override <- resolve_mask(df, feature_not_numeric)

  # Detection Logic
  if (is.null(cat_override)) {
    is_cat <- vapply(seq_len(n_cols), function(i) {
      col <- df[[i]]
      if (is.numeric(col) || is.integer(col)) return(FALSE)

      # Heuristic: if coercion creates *new* NAs, it's categorical
      ch  <- as.character(col)
      na0 <- is.na(ch)
      num <- suppressWarnings(as.numeric(ch))
      na1 <- is.na(num)
      any(na1 & !na0)
    }, logical(1))
  } else {
    is_cat <- as.logical(cat_override)
  }

  num_idx <- which(!is_cat)
  cat_idx <- which(is_cat)

  # Construct Numeric Matrix
  if (length(num_idx)) {
    # Efficient conversion
    mat_num <- tryCatch(
      as.matrix(df[, num_idx, drop = FALSE]),
      error = function(e) NULL
    )
    if (!is.numeric(mat_num)) {
      # Fallback coercion
      mat_num <- matrix(
        vapply(df[, num_idx, drop = FALSE], function(x) suppressWarnings(as.numeric(x)), numeric(n_rows)),
        nrow = n_rows
      )
      colnames(mat_num) <- colnames(df)[num_idx]
    }
  } else {
    mat_num <- matrix(numeric(0), nrow = n_rows, ncol = 0L)
  }

  # Construct Categorical Matrix
  if (length(cat_idx)) {
    mat_cat <- as.matrix(df[, cat_idx, drop = FALSE])
    mode(mat_cat) <- "character"
  } else {
    mat_cat <- matrix(character(0), nrow = n_rows, ncol = 0L)
  }

  list(
    numeric_data = mat_num,
    categoric_data = mat_cat,
    NumericData = mat_num,    # Compatibility alias
    CategoricData = mat_cat,  # Compatibility alias
    CategoricalData = mat_cat # Compatibility alias
  )
}

#' Split predictors into numeric and categorical matrices
#'
#' @description
#' Public wrapper that splits a dataset into numeric and categorical matrices
#' using the same logic as the internal helper used by tbc(). Prefer relying on
#' tbc()'s native handling; this is provided for convenience.
#'
#' @param data A matrix or data.frame.
#' @param feature_not_numeric Optional mask (logical/integer/character) indicating
#'   categorical columns.
#' @return A list with `numeric_data` and `categoric_data` matrices (and legacy
#'   compatibility aliases `NumericData`, `CategoricData`, `CategoricalData`).
#' @export
split_numeric_categoric <- function(data, feature_not_numeric = NULL) {
  .tbc_split_numeric_categoric(data, feature_not_numeric)
}

#' Split Predictors into Numeric and Categorical Matrices (Deprecated)
#'
#' @description
#' \strong{Deprecated}: Use internal logic or \code{tbc_fit}'s native handling instead.
#' This function remains for backward compatibility.
#'
#' @param Data Matrix or data.frame.
#' @param FeatureNotNumeric Categorical column selector.
#' @export
SplitDatasetInNumericAndCategoric <- function(Data, FeatureNotNumeric = NULL) {
  .Deprecated(".tbc_split_numeric_categoric")
  .tbc_split_numeric_categoric(Data, FeatureNotNumeric)
}

#' Standardize Column Names (Deprecated)
#' @keywords internal
standardize_colnames4plotting <- function(Header) {
  # Deprecated logic kept for specific legacy plotters if any
  Header <- gsub("\\.a", "\\.+", Header)
  Header <- gsub("\\.n", "\\.-", Header)
  gsub("\\.\\.", " ", Header)
}

#' Apply Categorical Mask to Data
#' @keywords internal
.tbc_apply_categoric_mask <- function(x, categoric_idx) {
  df <- as.data.frame(x, stringsAsFactors = FALSE, optional = TRUE)
  if (is.null(names(df))) names(df) <- paste0("V", seq_len(ncol(df)))

  if (is.null(categoric_idx)) return(df)

  mask <- categoric_idx
  if (is.logical(mask)) {
    mask <- which(mask)
  } else {
    mask <- as.integer(mask)
  }

  mask <- mask[mask >= 1L & mask <= ncol(df)]
  for (j in mask) df[[j]] <- as.character(df[[j]])

  df
}

#' Validate Control Parameters
#' @keywords internal
.tbc_validate_control <- function(control, context = "tbc()") {
  if (!is.list(control)) {
    stop(context, ": 'control' must be a list.", call. = FALSE)
  }

  check_num <- function(val, name, int = FALSE, lower = -Inf, allow_inf = FALSE) {
    if (is.null(val)) return()
    if (!is.numeric(val) || length(val) != 1L) stop(context, ": ", name, " must be scalar numeric.", call. = FALSE)
    if (!allow_inf && !is.finite(val)) stop(context, ": ", name, " must be finite.", call. = FALSE)
    if (val < lower) stop(context, ": ", name, " must be >= ", lower, ".", call. = FALSE)
    if (int && abs(val - round(val)) > sqrt(.Machine$double.eps)) stop(context, ": ", name, " must be integer.", call. = FALSE)
  }

  check_bool <- function(val, name) {
    if (is.null(val)) return()
    if (!is.logical(val) || length(val) != 1L || is.na(val)) stop(context, ": ", name, " must be TRUE/FALSE.", call. = FALSE)
  }

  check_num(control$min_split, "control$min_split", int = TRUE, lower = 1)
  check_num(control$max_depth, "control$max_depth", int = TRUE, lower = 1)
  check_num(control$max_cat_levels, "control$max_cat_levels", int = TRUE, lower = 1)

  check_bool(control$bivariate, "control$bivariate")
  check_bool(control$shuffle_features, "control$shuffle_features")
  check_bool(control$pruning, "control$pruning")
  check_bool(control$debug, "control$debug")

  if (!is.null(control$seed)) check_num(control$seed, "control$seed", int = TRUE)

  if (!is.null(control$prior_probs)) {
    if (!is.numeric(control$prior_probs) || any(control$prior_probs < 0)) {
      stop(context, ": prior_probs must be non-negative numeric.", call. = FALSE)
    }
  }

  invisible(control)
}

#' Validate Training Data Dimensions
#' @keywords internal
.tbc_validate_training_data <- function(x, y, context = "tbc()") {
  if (is.null(x)) stop(context, ": 'x' cannot be NULL.", call. = FALSE)

  n <- nrow(x)
  if (n == 0L) stop(context, ": 'x' has 0 rows.", call. = FALSE)
  if (ncol(x) == 0L) stop(context, ": 'x' has 0 columns.", call. = FALSE)

  if (is.null(y)) stop(context, ": 'y' cannot be NULL.", call. = FALSE)
  if (length(y) != n) stop(context, ": length(y) must match nrow(x).", call. = FALSE)
  if (all(is.na(y))) stop(context, ": 'y' contains only NAs.", call. = FALSE)

  invisible(NULL)
}