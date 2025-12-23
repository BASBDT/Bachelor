#' Compute optimal pruning sequence for a TreeList
#'
#' Adds the optimal cost-complexity pruning sequence (CART-style) to a
#' TreeList produced by the internal learning functions.
#'
#' @param tree_list A list containing the tree structure.
#' @return The modified \code{tree_list} with pruning fields added.
#' @keywords internal
.tbc_prune_sequence <- function(tree_list) {

  if (missing(tree_list) || !is.list(tree_list)) {
    stop("tbc: internal error: 'tree_list' must be a list.", call. = FALSE)
  }

  required <- c("NodeNo", "Parent_NodeNo", "Children_NodeNo", "Is_Leaf", "Risk")
  if (!all(required %in% names(tree_list))) {
    stop("tbc: internal error: 'tree_list' is missing required components.", call. = FALSE)
  }

  n_nodes    <- length(tree_list$NodeNo)
  parent_ids <- as.integer(tree_list$Parent_NodeNo)
  child_ids  <- as.matrix(tree_list$Children_NodeNo)
  is_leaf    <- as.logical(tree_list$Is_Leaf)
  node_risk  <- as.numeric(tree_list$Risk)

  branch_risk     <- rep(0.0, n_nodes)
  n_leaves_branch <- rep(0L, n_nodes)

  leaf_indices <- which(is_leaf)
  branch_risk[leaf_indices]     <- node_risk[leaf_indices]
  n_leaves_branch[leaf_indices] <- 1L

  for (i in seq(n_nodes, 1L, by = -1L)) {
    if (!is_leaf[i]) {
      left  <- child_ids[i, 1L]
      right <- child_ids[i, 2L]

      if (!is.na(left) && left > 0L) {
        branch_risk[i]     <- branch_risk[i] + branch_risk[left]
        n_leaves_branch[i] <- n_leaves_branch[i] + n_leaves_branch[left]
      }
      if (!is.na(right) && right > 0L) {
        branch_risk[i]     <- branch_risk[i] + branch_risk[right]
        n_leaves_branch[i] <- n_leaves_branch[i] + n_leaves_branch[right]
      }
    }
  }

  when_pruned <- integer(n_nodes)
  max_steps   <- max(1L, sum(!is_leaf)) + 10L
  all_alpha   <- numeric(max_steps)
  n_term      <- integer(max_steps)

  prune_step  <- 0L
  all_alpha[1L] <- 0.0
  n_term[1L]    <- n_leaves_branch[1L]

  branches <- which(!is_leaf)

  while (length(branches) > 0L) {
    denom <- pmax(1e-10, n_leaves_branch[branches] - 1.0)
    numer <- node_risk[branches] - branch_risk[branches]
    g_val <- numer / denom

    min_g    <- min(g_val)
    best_idx <- which(g_val <= min_g + 1e-12)
    to_prune <- branches[best_idx]

    best_alpha <- min_g
    if (prune_step > 0L && best_alpha < all_alpha[prune_step + 1L]) {
      best_alpha <- all_alpha[prune_step + 1L]
    }

    prune_step <- prune_step + 1L

    new_leaves <- to_prune
    when_pruned[new_leaves] <- prune_step

    if (length(new_leaves) > 0L) {
      for (node_idx in new_leaves) {
        diff_cost  <- node_risk[node_idx] - branch_risk[node_idx]
        diff_count <- n_leaves_branch[node_idx] - 1L

        curr <- node_idx
        while (curr != 0L) {
          parent <- parent_ids[curr]
          if (!is.na(parent) && parent > 0L) {
            n_leaves_branch[parent] <- n_leaves_branch[parent] - diff_count
            branch_risk[parent]     <- branch_risk[parent] + diff_cost
            curr <- parent
          } else {
            curr <- 0L
          }
        }
        n_leaves_branch[node_idx] <- 1L
        branch_risk[node_idx]     <- node_risk[node_idx]
      }
    }

    all_alpha[prune_step + 1L] <- best_alpha
    n_term[prune_step + 1L]    <- n_leaves_branch[1L]

    branches <- which(!is_leaf & n_leaves_branch > 1L)
  }

  final_len <- prune_step + 1L
  tree_list$prunelist  <- when_pruned
  tree_list$alpha      <- all_alpha[seq_len(final_len)]
  tree_list$ntermnodes <- n_term[seq_len(final_len)]

  tree_list
}

#' @title Legacy: Get Optimal Prune Sequence
#' @description
#' Deprecated alias for \code{.tbc_prune_sequence}.
#'
#' @inheritParams .tbc_prune_sequence
#' @export
GetOptimalPruneSequence <- function(Tree) {
  .Deprecated(".tbc_prune_sequence")
  .tbc_prune_sequence(tree_list = Tree)
}