#' Stratified k-fold cross-validation for TBC
#'
#' @description
#' Runs stratified K-fold (optionally repeated) CV using the TBC S3 interface.
#' Computes classification metrics without external dependencies.
#'
#' @param X data.frame or matrix of features.
#' @param Cls factor (preferred) or vector coercible to factor.
#' @param K integer, number of folds (default 5).
#' @param repeats integer, number of CV repetitions (default 1).
#' @param control a control list from \code{\link{tbc_control}}. Passed to each fit.
#' @param seed optional integer for reproducibility of fold assignment.
#' @param return_predictions logical; if TRUE, include per-fold predictions.
#' @return A list with:
#' \itemize{
#' \item \code{folds}: data.frame with per-(repeat,fold) metrics (acc, bacc, macro_f1, mcc).
#' \item \code{summary}: data.frame with mean and sd for each metric across all runs.
#' \item \code{predictions} (optional): list of data.frames with columns
#'       \code{replicate}, \code{fold}, \code{row}, \code{y_true}, \code{y_pred}.
#' }
#' @examples
#' \dontrun{
#' res <- tbc_cv(iris[,1:4], iris$Species,
#'               K = 3, control = tbc_control(pruning = FALSE, splitmin = 2))
#' res$summary
#' }
#' @export
tbc_cv <- function(X, Cls, K = 5L, repeats = 1L,
                   control = tbc_control(),
                   seed = NULL,
                   return_predictions = FALSE) {
  stopifnot(K >= 2L, repeats >= 1L)
  y <- as.factor(Cls)
  n <- length(y)
  if (!is.null(seed)) set.seed(seed)

  # stratified folds per repeat
  folds <- .tbc_stratified_folds(y, K, repeats)

  fold_rows <- integer(0)
  fold_rep  <- integer(0)
  fold_id   <- integer(0)
  fold_acc  <- numeric(0)
  fold_bacc <- numeric(0)
  fold_f1   <- numeric(0)
  fold_mcc  <- numeric(0)

  pred_list <- if (isTRUE(return_predictions)) list() else NULL
  Xdf <- if (is.data.frame(X)) X else as.data.frame(X, check.names = FALSE)

  for (r in seq_len(repeats)) {
    fold_assign <- folds[[r]]
    for (k in seq_len(K)) {
      test_idx  <- which(fold_assign == k)
      train_idx <- which(fold_assign != k)

      fit <- tbc(Xdf[train_idx, , drop = FALSE], Cls = y[train_idx], control = control)
      pred <- predict(fit, Xdf[test_idx, , drop = FALSE], type = "class")

      # metrics: .tbc_metric_vec(truth, pred) -> acc, bacc, macro_f1, mcc
      m <- .tbc_metric_vec(y[test_idx], pred)

      fold_rows <- c(fold_rows, length(test_idx))
      fold_rep  <- c(fold_rep,  r)
      fold_id   <- c(fold_id,   k)
      fold_acc  <- c(fold_acc,  m["acc"])
      fold_bacc <- c(fold_bacc, m["bacc"])
      fold_f1   <- c(fold_f1,   m["macro_f1"])
      fold_mcc  <- c(fold_mcc,  m["mcc"])

      if (!is.null(pred_list)) {
        pred_list[[length(pred_list) + 1L]] <- data.frame(
          replicate = r,
          fold      = k,
          row       = test_idx,
          y_true    = y[test_idx],
          y_pred    = pred,
          stringsAsFactors = FALSE
        )
      }
    }
  }

  folds_df <- data.frame(
    replicate = fold_rep,
    fold      = fold_id,
    n_test    = fold_rows,
    acc       = fold_acc,
    bacc      = fold_bacc,
    macro_f1  = fold_f1,
    mcc       = fold_mcc
  )

  summ <- data.frame(
    metric = c("acc", "bacc", "macro_f1", "mcc"),
    mean   = c(
      mean(fold_acc,  na.rm = TRUE),
      mean(fold_bacc, na.rm = TRUE),
      mean(fold_f1,   na.rm = TRUE),
      mean(fold_mcc,  na.rm = TRUE)
    ),
    sd     = c(
      stats::sd(fold_acc,  na.rm = TRUE),
      stats::sd(fold_bacc, na.rm = TRUE),
      stats::sd(fold_f1,   na.rm = TRUE),
      stats::sd(fold_mcc,  na.rm = TRUE)
    )
  )

  out <- list(folds = folds_df, summary = summ)
  if (!is.null(pred_list)) out$predictions <- pred_list
  out
}

#' Nested CV for simple hyperparameter selection
#'
#' @description
#' Performs nested CV: an outer K_o split for evaluation and an inner K_i split
#' to select hyperparameters from a small grid. Selection criterion defaults to
#' balanced accuracy.
#'
#' @param X,Cls as in \code{\link{tbc_cv}}.
#' @param outer_k integer, outer folds (default 5).
#' @param inner_k integer, inner folds (default 3).
#' @param param_grid data.frame with one column per control field you want to tune.
#'   Each row is a candidate config; values are assigned into \code{\link{tbc_control}}.
#' @param base_control control list providing defaults (see \code{\link{tbc_control}}); fields in \code{param_grid}
#'   override these per candidate.
#' @param seed optional base seed; each outer fold sets a derived seed.
#' @param score one of \code{"bacc"}, \code{"acc"}, \code{"macro_f1"}, \code{"mcc"}.
#' @return A list with:
#' \itemize{
#' \item \code{outer}: data.frame per outer fold with metrics and the selected index
#'       of \code{param_grid} used.
#' \item \code{selected_params}: list of control lists (length = outer_k).
#' \item \code{grid}: the \code{param_grid} you provided (for reference).
#' }
#' @examples
#' \dontrun{
#' grid <- expand.grid(
#'   splitmin = c(2, 10),
#'   bivariate = c(FALSE, TRUE),
#'   stringsAsFactors = FALSE
#' )
#' nc <- tbc_ncv(iris[,1:4], iris$Species,
#'               outer_k = 3, inner_k = 3,
#'               param_grid = grid,
#'               base_control = tbc_control(pruning = FALSE, shuffle_features = FALSE),
#'               seed = 1)
#' nc$outer
#' }
#' @export
tbc_ncv <- function(X, Cls,
                    outer_k = 5L,
                    inner_k = 3L,
                    param_grid,
                    base_control = tbc_control(),
                    seed = NULL,
                    score = c("bacc", "acc", "macro_f1", "mcc")) {
  score <- match.arg(score)
  stopifnot(is.data.frame(param_grid), nrow(param_grid) >= 1L)

  y <- as.factor(Cls)
  Xdf <- if (is.data.frame(X)) X else as.data.frame(X, check.names = FALSE)
  n  <- nrow(Xdf)
  if (!is.null(seed)) set.seed(seed)

  # Outer stratified folds
  folds_outer <- .tbc_stratified_folds(y, outer_k, repeats = 1L)[[1L]]

  selected <- vector("list", length = outer_k)
  res_outer <- data.frame(
    fold      = integer(0),
    n_train   = integer(0),
    n_test    = integer(0),
    acc       = numeric(0),
    bacc      = numeric(0),
    macro_f1  = numeric(0),
    mcc       = numeric(0),
    grid_index = integer(0)
  )

  # helper: build a control from base + row of param_grid
  build_control <- function(row_idx) {
    ctrl <- base_control
    for (nm in names(param_grid)) {
      ctrl[[nm]] <- param_grid[[nm]][[row_idx]]
    }
    ctrl
  }

  for (k in seq_len(outer_k)) {
    test_idx  <- which(folds_outer == k)
    train_idx <- which(folds_outer != k)

    # Inner CV for selection
    inner_res <- data.frame(idx = seq_len(nrow(param_grid)), score = NA_real_)
    for (gi in seq_len(nrow(param_grid))) {
      ctrl_g <- build_control(gi)
      cv_g <- tbc_cv(
        Xdf[train_idx, , drop = FALSE],
        y[train_idx],
        K        = inner_k,
        repeats  = 1L,
        control  = ctrl_g,
        seed     = if (is.null(seed)) NULL else seed + 1000L + gi,
        return_predictions = FALSE
      )
      inner_res$score[gi] <- cv_g$summary$mean[cv_g$summary$metric == score]
    }
    best_idx <- which.max(inner_res$score)
    selected[[k]] <- build_control(best_idx)

    # Fit on full inner-train and evaluate on outer-test
    fit <- tbc(Xdf[train_idx, , drop = FALSE], Cls = y[train_idx], control = selected[[k]])
    pred <- predict(fit, Xdf[test_idx, , drop = FALSE], type = "class")
    mv   <- .tbc_metric_vec(y[test_idx], pred)

    res_outer <- rbind(
      res_outer,
      data.frame(
        fold      = k,
        n_train   = length(train_idx),
        n_test    = length(test_idx),
        acc       = mv["acc"],
        bacc      = mv["bacc"],
        macro_f1  = mv["macro_f1"],
        mcc       = mv["mcc"],
        grid_index = best_idx
      )
    )
  }

  list(
    outer = res_outer,
    selected_params = selected,
    grid = param_grid
  )
}

# --- helpers -----------------------------------------------------------------

#' @keywords internal
.tbc_stratified_folds <- function(y, K, repeats = 1L) {
  y <- as.factor(y)
  n <- length(y)
  out <- vector("list", length = repeats)
  for (r in seq_len(repeats)) {
    # stratify per class
    idx <- integer(n)
    for (lev in levels(y)) {
      ids <- which(y == lev)
      # randomize within class
      ids <- sample(ids, length(ids), replace = FALSE)
      # distribute across folds
      split <- rep(seq_len(K), length.out = length(ids))
      idx[ids] <- split
    }
    out[[r]] <- idx
  }
  out
}
