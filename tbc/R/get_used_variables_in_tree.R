#' Get the variables used in a TBC tree
#'
#' @description
#' Returns the set of variables used by any split in the tree, with both indices
#' and names. Works for numeric, categorical, and derived (bivariate) splits.
#'
#' @param TreeList A TBC tree list (e.g., \code{fit$tree}).
#' @param include_derived Logical; if \code{TRUE}, include the pretty names of
#'   derived splits in \code{$derived_names}.
#' @return A list with:
#' \itemize{
#' \item \code{numeric_indices}: integer vector of numeric column indices used.
#' \item \code{categorical_indices}: integer vector of categorical column indices used
#'   (1-based within the \emph{categorical} block).
#' \item \code{base_names}: unique base (original) variable names used by any split.
#' \item \code{derived_names}: unique pretty names for derived splits (may be empty).
#' \item \code{names_all}: union of \code{base_names} and \code{derived_names}.
#' }
#' @examples
#' \dontrun{
#' fit <- tbc(iris[,1:4], Cls = iris$Species,
#'            control = tbc_control(bivariate = TRUE, pruning = FALSE, seed = 1))
#' vars <- GetUsedVariablesInTree(fit$tree)
#' vars$names_all
#' }
#' @keywords internal
GetUsedVariablesInTree <- function(TreeList, include_derived = TRUE) {
  get_used_variables_in_tree(tree_list = TreeList, include_derived = include_derived)
}

#' Variables used in a TBC tree (indices and names)
#'
#' @description
#' Returns the set of variables used by any split in the tree, with both indices
#' and names. Works for numeric, categorical, and derived (bivariate) splits.
#'
#' @param tree_list A TBC tree list (e.g., \code{fit$tree}).
#' @param include_derived Logical; if \code{TRUE}, include the pretty names of
#'   derived splits in \code{$derived_names}.
#' @return A list with:
#' \itemize{
#' \item \code{numeric_indices}: integer vector of numeric column indices used.
#' \item \code{categorical_indices}: integer vector of categorical column indices used
#'   (1-based within the \emph{categorical} block).
#' \item \code{base_names}: unique base (original) variable names used by any split.
#' \item \code{derived_names}: unique pretty names for derived splits (may be empty).
#' \item \code{names_all}: union of \code{base_names} and \code{derived_names}.
#' }
#' @examples
#' \dontrun{
#' fit <- tbc(iris[, 1:4], Cls = iris$Species,
#'            control = tbc_control(bivariate = TRUE, pruning = FALSE, seed = 1))
#' vars <- get_used_variables_in_tree(fit$tree)
#' vars$names_all
#' }
#' @export
get_used_variables_in_tree <- function(tree_list, include_derived = TRUE) {
  stopifnot(is.list(tree_list))
  dname  <- tree_list$DecisionPlaneName
  kids   <- tree_list$Children_NodeNo
  NDIM   <- tree_list$TreeInfo$NDIM %||% NA_integer_
  CDIM   <- tree_list$TreeInfo$CDIM %||% 0L
  Names  <- tree_list$VarNames %||% character()

  is_leaf <- (dname == 0L) | (kids[, 1] == 0L & kids[, 2] == 0L)
  split_nodes <- which(!is_leaf)

  num_idx <- integer()
  cat_idx <- integer()
  base_names <- character()
  derived_names <- character()

  # Helper: derived label consistent with .tbc_treelist_to_rules
  pretty_derived <- function(i, j, a, spec = NULL) {
    # Prefer an explicit pretty name if provided by the builder
    if (!is.null(spec) && !is.null(spec$name) && nzchar(spec$name)) {
      return(spec$name)
    }
    # Else mirror the exact expression style used in rules
    ni <- Names[i] %||% paste0("C", i)
    nj <- Names[j] %||% paste0("C", j)
    coef <- format(a, trim = TRUE)
    sprintf("(%s + %s*%s)", ni, coef, nj)
  }

  for (n in split_nodes) {
    idx <- dname[n]
    # try to get a human-readable name from the DecisionPlaneName names() attribute
    named_here <- NULL
    nm_attr <- names(dname)
    if (!is.null(nm_attr) && length(nm_attr) >= n) {
      candidate <- nm_attr[[n]]
      if (!is.null(candidate) && nzchar(candidate) && candidate != "Leafnode") {
        named_here <- candidate
      }
    }
    # Derived?
    if (!is.null(tree_list$is_derived) && isTRUE(tree_list$is_derived[[n]])) {
      sp <- tree_list$DerivedSpec[[n]]
      if (!is.null(sp$i) && !is.null(sp$j)) {
        num_idx <- c(num_idx, as.integer(sp$i), as.integer(sp$j))
        if (length(Names)) {
          base_names <- c(base_names, Names[sp$i], Names[sp$j])
          if (isTRUE(include_derived)) {
            derived_names <- c(derived_names, pretty_derived(sp$i, sp$j, sp$a, sp))
          }
        }
      }
      next
    }
    # Categorical?
    if (!is.na(idx) && idx < 0L && is.finite(NDIM)) {
      real_idx <- abs(idx) + NDIM
      cat_block_idx <- real_idx - NDIM
      cat_idx <- c(cat_idx, cat_block_idx)
      if (length(Names) >= real_idx) {
        base_names <- c(base_names, Names[real_idx])
      } else if (!is.null(named_here)) {
        base_names <- c(base_names, named_here)
      }
      next
    }
    # Numeric?
    if (!is.na(idx) && idx > 0L) {
      num_idx <- c(num_idx, idx)
      if (length(Names) >= idx) {
        base_names <- c(base_names, Names[idx])
      } else if (!is.null(named_here)) {
        base_names <- c(base_names, named_here)
      }
    }
  }

  num_idx <- sort(unique(num_idx))
  cat_idx <- sort(unique(cat_idx))
  base_names <- sort(unique(base_names))
  derived_names <- sort(unique(derived_names))
  names_all <- sort(unique(c(base_names, derived_names)))

  list(
    numeric_indices = num_idx,
    categorical_indices = cat_idx,
    base_names = base_names,
    derived_names = derived_names,
    names_all = names_all
  )
}
