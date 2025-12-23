#' Predict from a raw TreeList
#'
#' Low-level helper to obtain predictions directly from a \code{TreeList} structure,
#' bypassing the \code{tbc_model} S3 object. This function uses the same routing
#' logic as \code{\link{predict.tbc_model}} and is intended for use with the
#' raw output of \code{\link{tbc_fit}}.
#'
#' @param tree_list A list structure returned by \code{\link{tbc_fit}}.
#' @param x_num Numeric matrix of predictors.
#' @param x_cat Categorical matrix of predictors.
#' @param classes Optional character vector of class labels. If \code{NULL},
#'   labels are inferred from \code{tree_list$ClassProbability}.
#' @param type Character string: "class" (default) or "prob".
#'
#' @return A vector of classes (if \code{type = "class"}) or a matrix of
#'   probabilities (if \code{type = "prob"}).
#'
#' @export
tbc_predict_raw <- function(tree_list,
                            x_num = NULL,
                            x_cat = NULL,
                            classes = NULL,
                            type = c("class", "prob")) {

  # --- 1) Validation ---------------------------------------------------------

  if (!is.list(tree_list)) {
    stop("tbc_predict_raw: 'tree_list' must be a list.", call. = FALSE)
  }

  type <- match.arg(type)

  probs_all <- tree_list$ClassProbability
  if (is.null(probs_all)) {
    stop("tbc_predict_raw: 'tree_list' is missing $ClassProbability.", call. = FALSE)
  }

  n_classes <- ncol(probs_all)
  cls_names <- colnames(probs_all)

  # Resolve Class Labels
  if (is.null(classes)) {
    if (!is.null(cls_names)) {
      classes <- cls_names
    } else {
      classes <- as.character(seq_len(n_classes))
    }
  } else {
    classes <- as.character(classes)
    if (length(classes) != n_classes) {
      stop(sprintf("tbc_predict_raw: length(classes) is %d, but tree has %d classes.",
                   length(classes), n_classes), call. = FALSE)
    }
  }

  # --- 2) Input Coercion -----------------------------------------------------

  # Ensure Numeric Data is Matrix (Double)
  if (!is.null(x_num)) {
    x_num <- as.matrix(x_num)
    storage.mode(x_num) <- "double"
  }

  # Ensure Categorical Data is Matrix (Character)
  if (!is.null(x_cat)) {
    if (is.data.frame(x_cat)) {
      # Safe coercion of data.frame factors/logicals to character
      x_cat <- as.matrix(as.data.frame(lapply(x_cat, as.character), stringsAsFactors = FALSE))
    } else {
      x_cat <- as.matrix(x_cat)
    }
    storage.mode(x_cat) <- "character"
  }

  # Validate Dimensions
  n_dim_num <- if (!is.null(x_num)) ncol(x_num) else 0L
  n_dim_cat <- if (!is.null(x_cat)) ncol(x_cat) else 0L

  if (n_dim_num + n_dim_cat == 0L) {
    stop("tbc_predict_raw: At least one of 'x_num' or 'x_cat' must be provided.", call. = FALSE)
  }

  # Align Row Counts (Handle NULLs by creating empty matrices)
  n_obs <- if (!is.null(x_num)) nrow(x_num) else nrow(x_cat)

  if (is.null(x_num)) {
    x_num <- matrix(0.0, nrow = n_obs, ncol = 0L)
  }
  if (is.null(x_cat)) {
    x_cat <- matrix("", nrow = n_obs, ncol = 0L)
  }

  if (nrow(x_num) != n_obs || nrow(x_cat) != n_obs) {
    stop("tbc_predict_raw: Row counts of 'x_num' and 'x_cat' must match.", call. = FALSE)
  }

  # --- 3) Prediction ---------------------------------------------------------

  # Delegate to the central routing logic (from tbc.R)
  node_ids <- .tbc_assign_nodes(tree_list, x_num, x_cat)

  if (length(node_ids) == 0L) {
    stop("tbc_predict_raw: No observations to predict.", call. = FALSE)
  }

  # Validate that observations ended up in valid leaves
  is_leaf <- tree_list$Is_Leaf
  if (!is.null(is_leaf)) {
    # Check if any assigned node is not marked as a leaf (sanity check)
    # Note: .tbc_assign_nodes should theoretically only return leaves,
    # but strictly checking ensures robustness against bad trees.
    if (any(!is_leaf[node_ids])) {
      stop("tbc_predict_raw: Some observations were assigned to non-leaf nodes.", call. = FALSE)
    }
  }

  # Extract Probabilities
  # We index directly into the probability matrix using node IDs
  probs_out <- probs_all[node_ids, , drop = FALSE]
  colnames(probs_out) <- classes

  if (type == "prob") {
    return(probs_out)
  }

  # Extract Classes
  idx_max <- max.col(probs_out, ties.method = "first")
  factor(classes[idx_max], levels = classes)
}