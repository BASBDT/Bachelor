#' Build Human-Readable Rules from TreeList
#'
#' Converts the internal TreeList structure into human-readable R expressions
#' for each leaf node. Handles numeric, categorical, and bivariate splits.
#'
#' @param tree_list List returned by \code{tbc_fit}.
#' @param include_na_routing Logical. If TRUE, explicit NA checks are embedded in rules.
#' @param derived_name Character. "pretty" (default) or "expr" for bivariate splits.
#'
#' @return A list with:
#' \item{RuleSet_R}{Named character vector of rules.}
#' \item{NamesInRuleSet}{Vector of feature names used.}
#' \item{LeafNo}{Integer vector of leaf node IDs.}
#' \item{ClassPerNode}{Integer vector of class predictions per leaf.}
#'
#' @keywords internal
.tbc_treelist_to_rules <- function(tree_list,
                                   include_na_routing = FALSE,
                                   derived_name = c("pretty", "expr")) {

  derived_name <- match.arg(derived_name)

  # --- 1) Validation & Setup -------------------------------------------------
  if (!is.list(tree_list)) stop("Invalid 'tree_list': must be a list.")
  if (is.null(tree_list$NodeNo) || is.null(tree_list$Children_NodeNo)) {
    stop("Invalid 'tree_list': missing NodeNo or Children_NodeNo.")
  }

  kids        <- tree_list$Children_NodeNo
  plane_names <- tree_list$DecisionPlaneName
  boundaries  <- tree_list$DecisionBoundary

  # Identify leaves: explicit 0-plane or 0-children
  is_leaf <- (plane_names == 0L) | (kids[, 1] == 0L & kids[, 2] == 0L)

  n_dim   <- if (!is.null(tree_list$TreeInfo$NDIM)) tree_list$TreeInfo$NDIM else NA_integer_
  v_names <- tree_list$VarNames %||% names(plane_names)
  na_left <- tree_list$NA_to_left %||% rep(TRUE, length(plane_names))

  # Accessors for Bivariate Splits
  has_derived <- !is.null(tree_list$DerivedSpec) && !is.null(tree_list$is_derived)

  get_derived <- function(node) {
    if (!has_derived) return(NULL)
    if (!isTRUE(tree_list$is_derived[[node]])) return(NULL)
    tree_list$DerivedSpec[[node]]
  }

  # Helper: Quote categorical levels for R syntax
  squote <- function(xs) {
    if (length(xs) == 0L) return(character())
    xs <- gsub("'", "\\\\'", xs, fixed = TRUE)
    paste0("'", xs, "'")
  }

  # --- 2) Variable Resolution ------------------------------------------------

  # Returns list(name = "display_name", expr = "R_code")
  resolve_var <- function(node) {
    # Case A: Bivariate (Derived)
    spec <- get_derived(node)
    if (!is.null(spec)) {
      i <- spec$i
      j <- spec$j
      a <- spec$a

      left <- if (!is.null(v_names) && length(v_names) >= max(i, j)) {
        list(ni = v_names[i], nj = v_names[j])
      } else {
        list(ni = paste0("C", i), nj = paste0("C", j))
      }

      expr <- sprintf("(%s + %s*%s)", left$ni, format(a, trim = TRUE), left$nj)

      shown <- if (derived_name == "pretty" && !is.null(spec$name) && nzchar(spec$name)) {
        spec$name
      } else {
        expr
      }
      return(list(name = shown, expr = expr))
    }

    # Case B: Standard (Numeric or Categorical)
    nm <- names(plane_names)[node]
    if (!is.null(nm) && nzchar(nm)) {
      return(list(name = nm, expr = nm))
    }

    # Fallback to index lookup
    idx <- as.integer(plane_names[node])
    if (!is.na(idx) && idx != 0L) {
      # Categorical index is negative -> map to absolute + n_dim offset
      real_idx <- if (idx < 0L && is.finite(n_dim)) abs(idx) + n_dim else idx

      nm2 <- if (!is.null(v_names) && length(v_names) >= real_idx) {
        v_names[real_idx]
      } else {
        paste0("V", real_idx)
      }
      return(list(name = nm2, expr = nm2))
    }

    list(name = "Leafnode", expr = "Leafnode")
  }

  # --- 3) Condition Builder --------------------------------------------------

  # Returns list(left = "str", right = "str", names_used = c(...))
  build_cond <- function(node) {
    v   <- resolve_var(node)
    bd  <- boundaries[[node]]
    nal <- isTRUE(na_left[node])
    rte <- isTRUE(include_na_routing)

    # --- Numeric Split ---
    if (is.numeric(bd) && length(bd) == 1L && is.finite(bd)) {
      cut <- format(bd, digits = 8, trim = TRUE)

      if (rte && nal) {
        # NA goes Left
        l_str <- sprintf("is.na(%s) | %s < %s", v$expr, v$expr, cut)
        r_str <- sprintf("!is.na(%s) & %s >= %s", v$expr, v$expr, cut)
      } else if (rte && !nal) {
        # NA goes Right
        l_str <- sprintf("!is.na(%s) & %s < %s", v$expr, v$expr, cut)
        r_str <- sprintf("is.na(%s) | %s >= %s", v$expr, v$expr, cut)
      } else {
        # No explicit NA routing in rule
        l_str <- sprintf("%s < %s", v$expr, cut)
        r_str <- sprintf("%s >= %s", v$expr, cut)
      }
      return(list(left = l_str, right = r_str, names_used = v$name))
    }

    # --- Categorical Split ---
    if (is.list(bd)) {
      A <- bd[[1]]
      B <- bd[[2]]
    } else {
      A <- bd
      B <- NULL
    }

    # Format sets
    fmt_set <- function(s) {
      s <- as.character(s[!is.na(s)])
      paste(squote(s), collapse = ", ")
    }

    A_str <- fmt_set(A)
    B_str <- if (!is.null(B)) fmt_set(B) else NULL

    # Logic construction
    in_A <- sprintf("%s %%in%% c(%s)", v$expr, A_str)

    # If B is known, use explicit membership, else use negation of A
    in_B <- if (!is.null(B_str)) {
      sprintf("%s %%in%% c(%s)", v$expr, B_str)
    } else {
      sprintf("!(%s %%in%% c(%s))", v$expr, A_str)
    }

    if (rte && nal) {
      l_str <- sprintf("is.na(%s) | %s", v$expr, in_A)
      r_str <- sprintf("!is.na(%s) & %s", v$expr, in_B)
    } else if (rte && !nal) {
      l_str <- sprintf("!is.na(%s) & %s", v$expr, in_A)
      r_str <- sprintf("is.na(%s) | %s", v$expr, in_B)
    } else {
      l_str <- in_A
      r_str <- in_B
    }

    return(list(left = l_str, right = r_str, names_used = v$name))
  }

  # --- 4) Tree Traversal (DFS) -----------------------------------------------

  rule_list  <- list()
  names_used <- character()
  classes    <- tree_list$ClassOfObjects %||% rep(NA_integer_, length(plane_names))

  dfs <- function(node, conds) {
    if (is_leaf[node]) {
      # Reached leaf: store rule
      leaf_expr <- if (length(conds) > 0) {
        paste0("which(", paste(conds, collapse = " & "), ")")
      } else {
        "which(rep(TRUE, length.out = nrow(X)))"
      }

      rule_list[[length(rule_list) + 1L]] <<- list(
        leaf  = node,
        class = classes[node],
        rule  = leaf_expr
      )
      return(invisible())
    }

    # Internal Node: recurse
    branch <- build_cond(node)
    names_used <<- c(names_used, branch$names_used)

    child_l <- kids[node, 1]
    child_r <- kids[node, 2]

    if (child_l > 0L) dfs(child_l, c(conds, branch$left))
    if (child_r > 0L) dfs(child_r, c(conds, branch$right))
  }

  # Start DFS
  dfs(1L, character())

  # --- 5) Output Construction ------------------------------------------------

  if (length(rule_list) == 0L) {
    return(list(
      RuleSet_R      = character(),
      NamesInRuleSet = character(),
      LeafNo         = integer(),
      ClassPerNode   = integer(),
      RuleSet_Prolog = NULL
    ))
  }

  # Format output vectors
  n_rules <- length(rule_list)
  vec_rules <- character(n_rules)
  vec_names <- character(n_rules)
  vec_leafs <- integer(n_rules)
  vec_class <- integer(n_rules)

  for (i in seq_len(n_rules)) {
    r <- rule_list[[i]]
    nm <- paste0("ind_C", r$class, "_", r$leaf)

    vec_rules[i] <- r$rule
    vec_names[i] <- nm
    vec_leafs[i] <- r$leaf
    vec_class[i] <- r$class
  }
  names(vec_rules) <- vec_names

  list(
    RuleSet_R      = vec_rules,
    NamesInRuleSet = unique(names_used[nzchar(names_used)]),
    LeafNo         = vec_leafs,
    ClassPerNode   = vec_class,
    RuleSet_Prolog = NULL
  )
}