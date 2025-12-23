# Internal: prune whole subtrees under the given branch nodes, compact, and remap indices
PruneTreeOnBranches <- function(Tree, Branches = NULL) {
  # no-op if nothing to prune
  if (is.null(Branches) || !length(Branches)) {
    # message("PruneTreeOnBranches: There are no branches to prune")
    return(Tree)
  }

  N <- length(Tree$NodeNo)

  # collect all descendants of the branch nodes
  parents <- as.integer(Branches)
  tokeep  <- rep(TRUE, N)
  kids    <- integer(0)

  repeat {
    newkids <- as.vector(Tree$Children_NodeNo[parents, , drop = FALSE])
    newkids <- newkids[newkids > 0L & !(newkids %in% kids)]
    if (!length(newkids)) break
    kids              <- c(kids, newkids)
    tokeep[newkids]   <- FALSE
    parents           <- newkids
  }

  # convert selected branches to leaves (remove split & children)
  idx <- as.integer(Branches)
  idx <- idx[idx >= 1L & idx <= N]
  if (!length(idx)) return(Tree)

  Tree$DecisionPlaneName[idx]  <- 0L
  # DecisionBoundary is a list → assign safely
  Tree$DecisionBoundary[idx]   <- rep(list(0), length(idx))
  Tree$Children_NodeNo[idx, ]  <- 0L

  # compact the tree by dropping pruned nodes, and build a node id remap
  ntokeep  <- sum(tokeep)
  nodemap  <- integer(N); nodemap[tokeep] <- seq_len(ntokeep)

  # subset all per-node vectors/matrices/lists
  Tree$Parent_NodeNo        <- Tree$Parent_NodeNo[tokeep]
  Tree$ClassOfObjects       <- Tree$ClassOfObjects[tokeep]
  Tree$DecisionPlaneName    <- Tree$DecisionPlaneName[tokeep]
  Tree$DecisionBoundary     <- Tree$DecisionBoundary[tokeep]
  Tree$Children_NodeNo      <- Tree$Children_NodeNo[tokeep, , drop = FALSE]
  Tree$NodeProbability      <- Tree$NodeProbability[tokeep]
  Tree$resubstitution_error <- Tree$resubstitution_error[tokeep]
  Tree$Risk                 <- Tree$Risk[tokeep]
  Tree$NodeSize             <- Tree$NodeSize[tokeep]
  Tree$ClassProbability     <- Tree$ClassProbability[tokeep, , drop = FALSE]
  Tree$ClassObjectNo        <- Tree$ClassObjectNo[tokeep, , drop = FALSE]
  Tree$NotInClass           <- Tree$NotInClass[tokeep]

  # recompute leaf flags and NodeNo
  Tree$Is_Leaf <- (Tree$DecisionPlaneName == 0L)
  Tree$NodeNo  <- seq_len(ntokeep)

  # remap parent/children indices to compacted ids
  maskP <- Tree$Parent_NodeNo > 0L
  Tree$Parent_NodeNo[maskP] <- nodemap[Tree$Parent_NodeNo[maskP]]

  maskC <- Tree$Children_NodeNo > 0L
  Tree$Children_NodeNo[maskC] <- nodemap[Tree$Children_NodeNo[maskC]]

  Tree
}


removeBadSplits=function(Tree){
  # Remove splits that contribute nothing to the tree.
  #
  N          = length(Tree$NodeNo)
  Is_Leaf    = Tree$Is_Leaf#t((Tree$DecisionPlaneName==0))                      # no split variable implies leaf NodeNo
  isntpruned = matrix(TRUE, 1, N)
  doprune    = matrix(FALSE, 1, N)
  Risk       = t(Tree$Risk)
  #adjfactor = (1 - 100*eps(ClassOfObjects(Risk)))
  adjfactor  = (1 - 100 * .Machine$double.eps)
  
  
  while(TRUE){                                                                  # Work up from the bottom of the tree
    # Find "twigs" with two leaf Children_NodeNo
    leafs      = which(Is_Leaf & isntpruned)                                      # elementwise
    branches   = which((!Is_Leaf) & isntpruned)#elementwise
    if(length(branches)==0){
      break
    }
    #twig = branches[apply(Is_Leaf[Tree$Children_NodeNo[branches,]], 1, sum) == 2]
    branchtemp = Tree$Children_NodeNo[branches,,drop=FALSE]
    
    twig       = branches[apply(cbind(branchtemp[, 1,drop=FALSE] %in% leafs,
                                      branchtemp[, 2,drop=FALSE] %in% leafs),
                                1, sum) == 2]
    
    if(length(twig) == 0){
      break
    }
    # must have just the root NodeNo left
    
    
    # Find twigs to "unsplit" if the error of the twig is no larger
    # than the sum of the errors of the Children_NodeNo
    Rtwig = Risk[twig]
    kids = Tree$Children_NodeNo[twig, ,drop=FALSE]
    Rsplit = apply(cbind(Risk[kids[, 1,drop=FALSE]], Risk[kids[, 2,drop=FALSE]]), 1, sum)
 
    unsplit = Rsplit >= Rtwig * adjfactor
    if (any(unsplit)) {
      # Mark Children_NodeNo as pruned, and mark twig as now a leaf
      isntpruned[as.vector(kids[unsplit, ,drop=FALSE])] = 0
      twig                                              = twig[unsplit]         # only these to be marked on next 2 lines
      Is_Leaf[twig]                                     = TRUE
      doprune[twig]                                     = TRUE
    }else{
      break
    }
  } #end while
  Tree$Is_Leaf = Is_Leaf
  # Remove splits that are useless
  Tree$UselessBranches = which(doprune)
  if(any(doprune)){
    Tree = PruneTreeOnBranches(Tree = Tree, Branches = which(doprune))
  }
  return(Tree)
}