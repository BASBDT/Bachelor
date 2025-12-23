SimpleAccuracy4Treelist <- function(Treelist) {
  tbc_treelist_accuracy(Treelist)
}

#' Training accuracy of a TreeList
#'
#' @description
#' Computes the simple training accuracy of a tree represented as a raw
#' TreeList. This is intended for diagnostics and may overestimate
#' generalization performance.
#'
#' @param tree_list A TreeList as returned by \code{tbc_fit()}.
#'
#' @return Numeric scalar between 0 and 1.
#' @export
tbc_treelist_accuracy <- function(tree_list) {
  if (is.matrix(tree_list$ClassObjectNo)) {
    number_per_node <- apply(tree_list$ClassObjectNo, 1L, sum)
  } else if (is.vector(tree_list$ClassObjectNo)) {
    number_per_node <- tree_list$ClassObjectNo
  } else {
    warning("tbc_treelist_accuracy: tree_list$ClassObjectNo inappropriate data structure found.")
    return(NA_real_)
  }

  leaf <- tree_list$Is_Leaf
  acc  <- 1 - sum(tree_list$NotInClass[leaf]) / sum(number_per_node[leaf])
  acc
}

#' Internal shim for legacy accuracy function name
#' @keywords internal
.tbc_calc_accuracy <- function(TreeList) {
  tbc_treelist_accuracy(tree_list = TreeList)
}