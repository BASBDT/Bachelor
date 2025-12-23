#' Leaf-wise Summary and Rule Extraction
#'
#' Internal utilities for analyzing the leaves of a TBC model, including
#' classification statistics and rule string generation.
#'
#' @name leaf_analysis
#' @keywords internal
NULL

#' Structure classification summary for leaves
#'
#' @param tree_list A TBC TreeList.
#' @return Matrix of leaf statistics.
#' @keywords internal
.tbc_leaf_stats_matrix <- function(tree_list) {
  parents <- tree_list$Parent_NodeNo
  leaves  <- which(tree_list$Is_Leaf)

  if (length(leaves) == 0L) {
    return(matrix(numeric(0), ncol = 5, dimnames = list(NULL, c("LeafNo", "Class", "NoCases", "NotInClass", "InClass"))))
  }

  # Calculate leaf cases (sum of rows in ClassObjectNo for leaf indices)
  # ClassObjectNo is N_nodes x N_classes
  cls_obj <- tree_list$ClassObjectNo

  if (is.null(cls_obj)) {
    stop(".tbc_leaf_stats_matrix: TreeList missing ClassObjectNo.", call. = FALSE)
  }

  leaf_cases <- rowSums(cls_obj[leaves, , drop = FALSE])
  leaf_err   <- tree_list$NotInClass[leaves]

  cbind(
    LeafNo     = leaves,
    Class      = tree_list$ClassOfObjects[leaves],
    NoCases    = leaf_cases,
    NotInClass = leaf_err,
    InClass    = leaf_cases - leaf_err
  )
}

#' Resolve TreeList from input
#' @keywords internal
.tbc_resolve_treelist <- function(object) {
  if (inherits(object, "tbc_model")) {
    tl <- object$tree
    if (is.null(tl)) stop("tbc_model has no $tree component.", call. = FALSE)
    return(tl)
  }
  if (is.list(object) && !is.null(object$NodeNo)) {
    return(object)
  }
  stop("Object must be either a tbc_model or a TreeList.", call. = FALSE)
}

#' Generate Leaf Summary Table
#'
#' Combines statistics and rule strings for every leaf.
#'
#' @param object A \code{tbc_model} or \code{TreeList}.
#' @param include_na_routing Logical. Include NA paths in rules?
#' @param derived_name Character. "pretty" or "expr".
#'
#' @return A data.frame with columns: leaf, class_id, class_label, n, in_class,
#'   not_in_class, majority_prob, rule.
#' @keywords internal
.tbc_leaf_summary <- function(object,
                              include_na_routing = TRUE,
                              derived_name = c("pretty", "expr")) {
  derived_name <- match.arg(derived_name)
  tl <- .tbc_resolve_treelist(object)

  # 1. Get Labels if available
  cls_labels <- NULL
  if (inherits(object, "tbc_model")) {
    cls_labels <- object$classes
  }

  # 2. Compute Statistics
  mat <- .tbc_leaf_stats_matrix(tl)
  df_stats <- as.data.frame(mat, stringsAsFactors = FALSE)

  # Map labels
  if (!is.null(cls_labels)) {
    # Check bounds
    max_cls <- max(df_stats$Class, na.rm = TRUE)
    if (length(cls_labels) < max_cls) {
      warning(".tbc_leaf_summary: Fewer labels than class IDs.")
      df_stats$class_label <- NA_character_
    } else {
      df_stats$class_label <- cls_labels[df_stats$Class]
    }
  } else {
    df_stats$class_label <- NA_character_
  }

  # Calc Purity
  df_stats$majority_prob <- ifelse(df_stats$NoCases > 0,
                                   df_stats$InClass / df_stats$NoCases,
                                   NA_real_)

  # 3. Generate Rules
  # Uses the refactored rule generator from Treelist2Rules.R
  rules_list <- .tbc_treelist_to_rules(
    tl,
    include_na_routing = include_na_routing,
    derived_name       = derived_name
  )

  # 4. Merge Rules with Stats
  # rules_list$LeafNo, rules_list$RuleSet_R
  # Note: .tbc_treelist_to_rules returns vectors in matching order

  df_rules <- data.frame(
    LeafNo = rules_list$LeafNo,
    rule   = unname(rules_list$RuleSet_R),
    stringsAsFactors = FALSE
  )

  merged <- merge(df_stats, df_rules, by = "LeafNo", all.x = TRUE, sort = FALSE)
  merged <- merged[order(merged$LeafNo), ]

  # Final formatting
  data.frame(
    leaf          = merged$LeafNo,
    class_id      = merged$Class,
    class_label   = merged$class_label,
    n             = merged$NoCases,
    in_class      = merged$InClass,
    not_in_class  = merged$NotInClass,
    majority_prob = merged$majority_prob,
    rule          = merged$rule,
    stringsAsFactors = FALSE
  )
}

#' Extract Top Rules for a Specific Class
#'
#' @param object A \code{tbc_model}.
#' @param class_target Target class (name or integer ID).
#' @param top_n Number of rules to return.
#' @param ... Passed to \code{.tbc_leaf_summary}.
#'
#' @return List with leaf table and rule vector.
#' @keywords internal
.tbc_extract_top_rules <- function(object, class_target, top_n = 5L, ...) {

  if (!inherits(object, "tbc_model")) {
    stop(".tbc_extract_top_rules: object must be a tbc_model.", call. = FALSE)
  }

  # Generate full summary
  leaf_tab <- .tbc_leaf_summary(object, ...)

  # Resolve Target ID
  if (is.character(class_target)) {
    idx <- match(class_target, object$classes)
    if (is.na(idx)) {
      stop("Class '", class_target, "' not found in model classes.", call. = FALSE)
    }
    target_id <- idx
  } else {
    target_id <- as.integer(class_target)
  }

  # Filter & Sort
  sub <- leaf_tab[leaf_tab$class_id == target_id, , drop = FALSE]

  if (nrow(sub) == 0L) {
    return(list(leaf_table = sub, top_rules = character(0)))
  }

  # Sort by Coverage (InClass desc), Purity (desc), Leaf ID (asc)
  ord <- order(-sub$in_class, -sub$majority_prob, sub$leaf)
  sub <- sub[ord, , drop = FALSE]

  k <- min(top_n, nrow(sub))
  top_rules <- sub$rule[seq_len(k)]

  list(
    leaf_table = sub,
    top_rules  = top_rules
  )
}

#' Analyze Rule Complexity
#'
#' @param rules_vec Character vector of rule strings.
#' @return List of complexity metrics.
#' @keywords internal
.tbc_count_rule_complexity <- function(rules_vec) {
  if (length(rules_vec) == 0L) {
    return(list(
      n_rules         = 0L,
      conds_per_rule  = integer(0),
      n_unique_conds  = 0L,
      n_total_conds   = 0L
    ))
  }

  # Helper to split "which(A & B & C)"
  split_conds <- function(r_str) {
    # Remove wrapper "which(...)"
    # Regex: remove everything up to first ( and last )
    inner <- sub("^which\\((.*)\\)$", "\\1", r_str)

    # Split by " & "
    parts <- unlist(strsplit(inner, " & ", fixed = TRUE))

    # Remove whitespace
    parts <- gsub(" ", "", parts)
    parts[nzchar(parts)]
  }

  conds_list <- lapply(rules_vec, split_conds)

  all_conds <- unlist(conds_list)

  list(
    n_rules         = length(rules_vec),
    conds_per_rule  = vapply(conds_list, length, integer(1)),
    n_unique_conds  = length(unique(all_conds)),
    n_total_conds   = length(all_conds)
  )
}