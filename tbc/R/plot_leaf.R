#' Plot Class Distribution for a Single Leaf
#'
#' Visualizes the training-time class frequencies within a specific leaf node.
#' If new \code{data} is provided, it counts how many of the new observations
#' fall into this leaf (using the tree's specific NA routing and split logic),
#' adding this count to the subtitle.
#'
#' @param object A \code{tbc_model} or a raw \code{TreeList}.
#' @param leaf_id Integer. The node ID of the leaf to plot.
#' @param data Optional data.frame or matrix. If provided, the number of rows
#'   routed to \code{leaf_id} is calculated and displayed.
#' @param include_na_routing Logical. Include NA routing logic in the displayed rule text?
#' @param derived_name Character. "pretty" or "expr" for bivariate split naming in the rule text.
#' @param ... Additional arguments passed to \code{\link[graphics]{barplot}}.
#'
#' @return Invisibly returns the bar midpoints (from \code{barplot}).
#' @export
tbc_plot_leaf <- function(object,
                          leaf_id,
                          data = NULL,
                          include_na_routing = TRUE,
                          derived_name = c("pretty", "expr"),
                          ...) {

  derived_name <- match.arg(derived_name)
  tl <- .tbc_resolve_treelist(object)
  leaf_id <- as.integer(leaf_id)

  # --- 1) Validate Leaf ID ---------------------------------------------------

  # Map external leaf_id to internal row index
  # Note: NodeNo usually equals 1:N, but we treat it safely as a lookup
  node_ids <- tl$NodeNo %||% seq_along(tl$DecisionPlaneName)
  idx <- match(leaf_id, node_ids)

  if (is.na(idx)) {
    stop("tbc_plot_leaf: Leaf ID ", leaf_id, " not found in tree.", call. = FALSE)
  }

  # Check if it is actually a leaf
  is_leaf <- if (!is.null(tl$Is_Leaf)) tl$Is_Leaf[idx] else (tl$DecisionPlaneName[idx] == 0)
  if (!is_leaf) {
    stop("tbc_plot_leaf: Node ", leaf_id, " is not a leaf node.", call. = FALSE)
  }

  # --- 2) Retrieve Leaf Statistics & Rule ------------------------------------

  # We use the internal summary helper to get the rule string and basic stats
  # efficiently without recalculating everything manually.
  leaf_df <- .tbc_leaf_summary(
    object,
    include_na_routing = include_na_routing,
    derived_name = derived_name
  )

  stats <- leaf_df[leaf_df$leaf == leaf_id, , drop = FALSE]

  if (nrow(stats) == 0L) {
    stop("tbc_plot_leaf: Could not retrieve summary statistics for leaf ", leaf_id, ".", call. = FALSE)
  }

  # --- 3) Get Training Counts ------------------------------------------------

  if (is.null(tl$ClassObjectNo)) {
    stop("tbc_plot_leaf: TreeList is missing 'ClassObjectNo' (training counts).", call. = FALSE)
  }

  # Extract counts for this node
  counts_train <- as.numeric(tl$ClassObjectNo[idx, , drop = TRUE])
  n_classes    <- length(counts_train)

  # Resolve Class Labels
  # Priority: 1. Column names of count matrix -> 2. Model classes -> 3. Generic
  cls_labels <- colnames(tl$ClassObjectNo)

  if (is.null(cls_labels) && inherits(object, "tbc_model")) {
    cls_labels <- object$classes
  }

  if (is.null(cls_labels) || length(cls_labels) != n_classes) {
    cls_labels <- paste0("Class ", seq_len(n_classes))
  }

  names(counts_train) <- cls_labels

  # --- 4) Handle New Data (Optional) -----------------------------------------

  n_new <- NA_integer_

  if (!is.null(data)) {
    # We must route the data to find out how many hit this leaf.
    # We reuse the robust logic from tbc_predict / tbc_predict_raw.

    if (inherits(object, "tbc_model")) {
      # Safe model-based preparation
      prep <- .tbc_prepare_newdata(object, data)
      x_num <- prep$x_num
      x_cat <- prep$x_cat
    } else {
      # Raw TreeList: best-effort coercion
      if (is.matrix(data) || is.data.frame(data)) {
        # Simple heuristic split based on data types
        split_res <- .tbc_split_data(as.data.frame(data, stringsAsFactors = FALSE))
        x_num <- split_res$x_num
        x_cat <- split_res$x_cat
      } else {
        warning("tbc_plot_leaf: 'data' provided but format unclear; skipping new data count.")
        x_num <- NULL; x_cat <- NULL
      }
    }

    if (!is.null(x_num) || !is.null(x_cat)) {
      assigned_nodes <- .tbc_assign_nodes(tl, x_num, x_cat)
      # We map the internal index 'idx' back to the external 'leaf_id'
      # But assigned_nodes returns the internal index (1..N).
      # Note: .tbc_assign_nodes returns *indices* into the node vectors,
      # which correspond to 'idx'. We need to check if assigned_nodes == idx.
      n_new <- sum(assigned_nodes == idx, na.rm = TRUE)
    }
  }

  # --- 5) Plotting -----------------------------------------------------------

  # Prepare Titles
  main_title <- sprintf("Leaf %d", leaf_id)

  # Majority info
  maj_cls   <- if (!is.na(stats$class_label)) stats$class_label else stats$class_id
  maj_prob  <- stats$majority_prob

  sub_title <- sprintf(
    "Train N = %d | Maj: %s (%.2f)",
    stats$n, maj_cls, maj_prob
  )

  if (!is.na(n_new)) {
    sub_title <- paste0(sub_title, sprintf(" | New Data N = %d", n_new))
  }

  rule_txt <- if (!is.na(stats$rule) && nzchar(stats$rule)) stats$rule else NULL

  # Save par settings
  op <- graphics::par(no.readonly = TRUE)
  on.exit(graphics::par(op), add = TRUE)

  # Adjust margins to fit rule text at bottom
  if (!is.null(rule_txt)) {
    graphics::par(mar = c(6, 4, 4, 2) + 0.1)
  }

  bp <- graphics::barplot(
    counts_train,
    main = main_title,
    ylab = "Training Count",
    xlab = "Class",
    las  = 1, # Horizontal axis labels
    ...
  )

  graphics::mtext(sub_title, side = 3, line = 0.5, cex = 0.8)

  if (!is.null(rule_txt)) {
    # Wrap text if too long
    wrapped_rule <- paste(strwrap(rule_txt, width = 60), collapse = "\n")
    graphics::mtext(wrapped_rule, side = 1, line = 4.5, cex = 0.7, adj = 0.5)
  }

  invisible(bp)
}