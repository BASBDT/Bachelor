# Internal: prune whole subtrees under the given branch nodes, compact, and remap indices
PruneTreeOnBranches <- function(Tree, Branches = NULL) {
  # no-op if nothing to prune
  if (is.null(Branches) || !length(Branches)) {
    return(Tree)
  }
  if (!is.list(Tree) || is.null(Tree$NodeNo)) {
    stop("PruneTreeOnBranches: `Tree` must be a TreeList with $NodeNo.", call. = FALSE)
  }

  N <- length(Tree$NodeNo)

  # sanitize branch indices (they are node indices in current representation)
  Branches <- as.integer(Branches)
  Branches <- Branches[is.finite(Branches) & Branches >= 1L & Branches <= N]
  if (!length(Branches)) {
    return(Tree)
  }

  # collect all descendants of the branch nodes
  parents <- Branches
  tokeep  <- rep(TRUE, N)
  kids    <- integer(0)

  repeat {
    child_mat <- Tree$Children_NodeNo[parents, , drop = FALSE]
    newkids   <- as.vector(child_mat)
    newkids   <- newkids[newkids > 0L]
    if (!length(newkids)) break

    # drop duplicates / already seen
    newkids <- newkids[!duplicated(newkids)]
    newkids <- newkids[!(newkids %in% kids)]
    if (!length(newkids)) break

    kids            <- c(kids, newkids)
    tokeep[newkids] <- FALSE
    parents         <- newkids
  }

  # convert selected branches themselves to leaves (remove split & children)
  idx <- Branches
  if (!length(idx)) return(Tree)

  Tree$DecisionPlaneName[idx] <- 0L
  if (!is.null(Tree$DecisionBoundary)) {
    Tree$DecisionBoundary[idx] <- rep(list(0), length(idx))
  }
  if (!is.null(Tree$Children_NodeNo)) {
    Tree$Children_NodeNo[idx, ] <- 0L
  }

  # compact the tree by dropping pruned descendant nodes, and build node id remap
  ntokeep <- sum(tokeep)
  nodemap <- integer(N)
  nodemap[tokeep] <- seq_len(ntokeep)

  subset_vec  <- function(x) if (!is.null(x)) x[tokeep] else NULL
  subset_mat  <- function(x) if (!is.null(x)) x[tokeep, , drop = FALSE] else NULL
  subset_list <- function(x) if (!is.null(x)) x[tokeep] else NULL

  # subset all per-node vectors/matrices/lists we maintain elsewhere
  Tree$Parent_NodeNo        <- subset_vec(Tree$Parent_NodeNo)
  Tree$ClassOfObjects       <- subset_vec(Tree$ClassOfObjects)
  Tree$DecisionPlaneName    <- subset_vec(Tree$DecisionPlaneName)
  Tree$DecisionBoundary     <- subset_list(Tree$DecisionBoundary)
  Tree$Children_NodeNo      <- subset_mat(Tree$Children_NodeNo)
  Tree$NodeProbability      <- subset_vec(Tree$NodeProbability)
  Tree$resubstitution_error <- subset_vec(Tree$resubstitution_error)
  Tree$Risk                 <- subset_vec(Tree$Risk)
  Tree$NodeSize             <- subset_vec(Tree$NodeSize)
  Tree$ClassProbability     <- subset_mat(Tree$ClassProbability)
  Tree$ClassObjectNo        <- subset_mat(Tree$ClassObjectNo)
  Tree$NotInClass           <- subset_vec(Tree$NotInClass)

  # newer per-node fields: NA routing and bivariate metadata
  if (!is.null(Tree$NA_to_left)) {
    Tree$NA_to_left <- Tree$NA_to_left[tokeep]
  }
  if (!is.null(Tree$is_derived)) {
    Tree$is_derived <- Tree$is_derived[tokeep]
  }
  if (!is.null(Tree$DerivedSpec)) {
    Tree$DerivedSpec <- Tree$DerivedSpec[tokeep]
  }
  if (!is.null(Tree$prunelist)) {
    Tree$prunelist <- Tree$prunelist[tokeep]
  }

  # recompute leaf flags and NodeNo
  Tree$Is_Leaf <- (Tree$DecisionPlaneName == 0L)
  Tree$NodeNo  <- seq_len(ntokeep)

  # remap parent/children indices to compacted ids
  if (!is.null(Tree$Parent_NodeNo)) {
    maskP <- Tree$Parent_NodeNo > 0L
    Tree$Parent_NodeNo[maskP] <- nodemap[Tree$Parent_NodeNo[maskP]]
  }

  if (!is.null(Tree$Children_NodeNo)) {
    maskC <- Tree$Children_NodeNo > 0L
    Tree$Children_NodeNo[maskC] <- nodemap[Tree$Children_NodeNo[maskC]]
  }

  Tree
}

#' Internal shim for legacy branch-pruning API
#' @keywords internal
.tbc_prune_branches <- function(tree_list, branches = NULL) {
  PruneTreeOnBranches(Tree = tree_list, Branches = branches)
}
