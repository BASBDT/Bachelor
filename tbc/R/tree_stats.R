# Tree size / depth utilities ----------------------------------------------

#' Depth per node for a TBC TreeList
#'
#' Computes the depth (distance from root) for each node in a TreeList.
#' Root depth is 1, its children have depth 2, and so on.
#'
#' @param TreeList A TreeList as returned by \code{LearnPlanarDecision()}.
#' @return An integer vector of length = number of nodes. Names are the
#'         node ids (\code{TreeList$NodeNo}) if available.
#' @keywords internal
TreeDepthPerNode <- function(TreeList) {
  if (is.null(TreeList$Children_NodeNo)) {
    stop("TreeDepthPerNode: TreeList has no 'Children_NodeNo' field.", call. = FALSE)
  }

  kids <- TreeList$Children_NodeNo
  if (!is.matrix(kids) || ncol(kids) != 2L) {
    stop("TreeDepthPerNode: 'Children_NodeNo' must be an n x 2 matrix.", call. = FALSE)
  }

  node_ids <- TreeList$NodeNo
  if (is.null(node_ids)) {
    node_ids <- seq_len(nrow(kids))
  }

  # Root = node that is never a child; fallback to first node if ambiguous
  all_children <- as.integer(kids)
  all_children <- all_children[is.finite(all_children) & all_children > 0L]
  roots <- setdiff(node_ids, all_children)
  root_id <- if (length(roots)) roots[1L] else node_ids[1L]

  idx_by_id <- seq_along(node_ids)
  names(idx_by_id) <- as.character(node_ids)

  root_idx <- idx_by_id[[as.character(root_id)]]
  if (is.na(root_idx)) {
    stop("TreeDepthPerNode: could not identify root node.", call. = FALSE)
  }

  depth <- rep(NA_integer_, length(node_ids))
  depth[root_idx] <- 1L

  queue <- root_idx
  while (length(queue)) {
    i <- queue[1L]
    queue <- queue[-1L]

    d <- depth[i]
    child_ids <- kids[i, ]
    for (cid in child_ids) {
      if (!is.finite(cid) || cid <= 0L) next
      j <- idx_by_id[[as.character(cid)]]
      if (is.na(j)) next
      if (is.na(depth[j]) || depth[j] > d + 1L) {
        depth[j] <- d + 1L
        queue <- c(queue, j)
      }
    }
  }

  names(depth) <- as.character(node_ids)
  depth
}

#' Tree size summary for TBC models or TreeLists
#'
#' @param object A \code{tbc_model} or a raw TreeList.
#' @return A list with elements \code{depth}, \code{n_nodes},
#'   \code{n_internal}, and \code{n_leaves}.
#' @export
tbc_tree_size <- function(object) {
  if (inherits(object, "tbc_model")) {
    TL <- object$tree
  } else if (is.list(object) && !is.null(object$Children_NodeNo)) {
    TL <- object
  } else {
    stop("tbc_tree_size: expected a tbc_model or a TreeList.", call. = FALSE)
  }

  kids <- TL$Children_NodeNo
  if (is.null(kids) || !is.matrix(kids) || ncol(kids) != 2L) {
    stop("tbc_tree_size: TreeList is missing a valid 'Children_NodeNo' matrix.", call. = FALSE)
  }

  depth_vec <- TreeDepthPerNode(TL)
  depth <- if (length(depth_vec)) max(depth_vec, na.rm = TRUE) else NA_integer_

  # Leaf detection: prefer explicit Is_Leaf, otherwise derive from children
  if (!is.null(TL$Is_Leaf)) {
    is_leaf <- as.integer(TL$Is_Leaf) != 0L
  } else {
    is_leaf <- (kids[, 1L] <= 0L & kids[, 2L] <= 0L)
  }

  n_nodes   <- length(depth_vec)
  n_leaves  <- sum(is_leaf)
  n_internal <- n_nodes - n_leaves

  list(
    depth      = depth,
    n_nodes    = n_nodes,
    n_internal = n_internal,
    n_leaves   = n_leaves
  )
}
