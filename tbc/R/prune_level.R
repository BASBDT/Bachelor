#' Prune Tree to a Given Cost-Complexity Level
#'
#' @description
#' Prunes a TreeList according to its optimal cost–complexity pruning sequence
#' (as computed by \code{\link{.tbc_prune_sequence}}). All nodes whose
#' pruning step \code{<= level} are removed. The resulting TreeList still
#' carries the pruning metadata (\code{prunelist}, \code{alpha},
#' \code{ntermnodes}) from \code{.tbc_prune_sequence()}.
#'
#' @param tree_list TreeList as returned from \code{LearnPlanarDecision()}.
#' @param level Integer (>= 1). Maximum pruning step to apply.
#'
#' @return A pruned TreeList with the same structure as \code{tree_list}, but with
#'         selected branches collapsed to leaves. The fields
#'         \code{$prunelist}, \code{$alpha}, and \code{$ntermnodes} are present
#'         as produced by \code{.tbc_prune_sequence()}.
#' @keywords internal
.tbc_prune_level <- function(tree_list, level = NULL) {

  # ---- 1. Validation --------------------------------------------------------
  if (is.null(level)) {
    stop("tbc: 'level' must be supplied (integer >= 1) in .tbc_prune_level().",
         call. = FALSE)
  }
  if (!is.numeric(level) || length(level) != 1L || is.na(level)) {
    stop("tbc: 'level' must be a numeric scalar in .tbc_prune_level().",
         call. = FALSE)
  }

  level <- as.integer(level)
  if (level < 1L) level <- 1L

  if (!is.list(tree_list) || is.null(tree_list$NodeNo)) {
    stop("tbc: 'tree_list' must be a TreeList (list with $NodeNo).",
         call. = FALSE)
  }

  # ---- 2. Pruning Sequence --------------------------------------------------
  # Compute / attach optimal prune sequence using internal helper
  # (Assumes GetOptimalPruneSequence is refactored to .tbc_prune_sequence)
  tree_seq <- .tbc_prune_sequence(tree_list = tree_list)
  when_pruned <- tree_seq$prunelist

  if (is.null(when_pruned) || length(when_pruned) == 0L) {
    warning("tbc: .tbc_prune_sequence() returned empty 'prunelist'; returning tree unchanged.",
            call. = FALSE)
    return(tree_seq)
  }

  # ---- 3. Branch Selection --------------------------------------------------
  # Select branches to prune at or before the requested level
  branches <- which(when_pruned > 0L & when_pruned <= level)

  if (length(branches) == 0L) {
    # nothing to prune at this level -> return tree with pruning metadata
    return(tree_seq)
  }

  # ---- 4. Apply Pruning -----------------------------------------------------
  # Apply structural pruning on these branches
  # (Assumes PruneTreeOnBranches is refactored to .tbc_prune_branches)
  .tbc_prune_branches(tree_list = tree_seq, branches = branches)
}

#' Internal shim for legacy prune-by-level API
#' @keywords internal
.tbc_prune_by_level <- function(Tree, Level = NULL) {
  .tbc_prune_level(tree_list = Tree, level = Level)
}

#' @title Legacy: Prune Tree by Level
#' @description
#' Deprecated alias for \code{.tbc_prune_level}.
#'
#' @param Tree Mapped to \code{tree_list}.
#' @param Level Mapped to \code{level}.
#' @export
PruneTreeByXLevel <- function(Tree, Level = NULL) {
  .Deprecated(".tbc_prune_level")
  .tbc_prune_level(tree_list = Tree, level = Level)
}