#' Row Indices per Leaf for a TBC Tree
#'
#' @description
#' Traverses the tree for each row in \code{data} and returns a list containing,
#' for each leaf, the vector of row indices that fall into that leaf.
#'
#' This function handles numeric, categorical, and derived (bivariate) splits
#' by delegating to the central routing logic.
#'
#' @param tree_list A TBC tree list (e.g., \code{fit$tree}).
#' @param data A data.frame or matrix of features.
#' @param var_names Optional feature names; defaults to \code{tree_list$VarNames}
#'   or \code{colnames(data)}.
#'
#' @return A named list: each element name is the leaf node ID, and the value
#'   is an integer vector of row indices reaching that leaf.
#'
#' @keywords internal
.tbc_leaf_indices <- function(tree_list, data, var_names = NULL) {

  # --- 1) Validation ---------------------------------------------------------

  if (!is.list(tree_list)) {
    stop(".tbc_leaf_indices: 'tree_list' must be a list.", call. = FALSE)
  }

  if (is.null(data)) {
    stop(".tbc_leaf_indices: 'data' is required.", call. = FALSE)
  }

  # Ensure data is a data.frame for name-based extraction
  df <- as.data.frame(data, stringsAsFactors = FALSE)

  # Resolve Feature Names
  if (is.null(var_names)) {
    var_names <- tree_list$VarNames %||% colnames(df)
  }

  # Resolve Dimensions from Tree Info
  # NDIM tells us the split point between numeric and categorical features
  # in the internal training order.
  n_dim_num <- tree_list$TreeInfo$NDIM %||% 0L

  # --- 2) Prepare Matrices (x_num / x_cat) -----------------------------------

  # The tree logic (and .tbc_assign_nodes) relies on split variables being
  # indices into specific numeric or categorical matrices. We must reconstruct
  # these matrices based on the tree's VarNames.

  n_total <- length(var_names)
  n_dim_cat <- n_total - n_dim_num

  x_num <- NULL
  x_cat <- NULL

  if (length(var_names) > 0) {

    # 2a. Numeric Features (Indices 1 .. n_dim_num)
    if (n_dim_num > 0) {
      names_num <- var_names[seq_len(n_dim_num)]

      # Extract by name if present, otherwise by index
      # We create a matrix of NA if columns are missing
      mat_num <- matrix(NA_real_, nrow = nrow(df), ncol = n_dim_num)
      colnames(mat_num) <- names_num

      found_num <- intersect(names_num, colnames(df))
      if (length(found_num) > 0) {
        # Suppress coercion warnings
        suppressWarnings({
          tmp <- as.matrix(df[, found_num, drop = FALSE])
          mode(tmp) <- "numeric"
        })
        mat_num[, found_num] <- tmp
      }
      x_num <- mat_num
    }

    # 2b. Categorical Features (Indices n_dim_num+1 .. n_total)
    if (n_dim_cat > 0) {
      names_cat <- var_names[seq(from = n_dim_num + 1, length.out = n_dim_cat)]

      mat_cat <- matrix(NA_character_, nrow = nrow(df), ncol = n_dim_cat)
      colnames(mat_cat) <- names_cat

      found_cat <- intersect(names_cat, colnames(df))
      if (length(found_cat) > 0) {
        # Convert to character matrix
        tmp <- as.matrix(df[, found_cat, drop = FALSE])
        mode(tmp) <- "character"
        mat_cat[, found_cat] <- tmp
      }
      x_cat <- mat_cat
    }
  } else {
    # Fallback: if no var_names, try to split purely by column type or order?
    # This is risky for a trained tree. We assume data provided matches training structure.
    # For robust routing, we need the structure. If var_names is empty, likely an edge case.
    # We pass NULLs to .tbc_assign_nodes which handles empty data gracefully (returns empty).
  }

  # --- 3) Route Observations -------------------------------------------------

  # Use the central routing helper
  # (Returns vector of node IDs, one per row)
  node_ids <- .tbc_assign_nodes(tree_list, x_num, x_cat)

  if (length(node_ids) == 0L) {
    return(list())
  }

  # --- 4) Group Indices by Leaf ----------------------------------------------

  # Identify actual leaf nodes in the tree structure
  # (Some rows might hypothetically be routed to 0 or non-leaves if tree is bad, though assign_nodes prevents this)

  # split() returns a list of indices
  row_indices <- seq_len(nrow(df))
  res <- split(row_indices, node_ids)

  # Filter out '0' node if it exists (dropped rows)
  if (!is.null(res$`0`)) {
    res$`0` <- NULL
  }

  res
}