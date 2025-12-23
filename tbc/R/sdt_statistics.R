#' Structural statistics for decision trees across engines
#'
#' Compute unified structural metrics for trees produced by several engines:
#' rpart, partykit/ctree/evtree, tree, C5.0, RWeka J48, LDATree (Treee), and
#' CORElearn (CoreModel with model = "tree"). Internally, where possible,
#' models are converted to `partykit::party` to use a common traversal; for
#' engines without a direct bridge we compute metrics from their native
#' structures.
#'
#' The returned list contains:
#' - `min_depth`: minimum leaf depth
#' - `avg_depth`: mean leaf depth
#' - `max_depth`: maximum leaf depth
#' - `avg_weighted_depth`: leaf-depth average weighted by number of training
#'   samples per leaf
#' - `decision_nodes`: number of internal (non-leaf) nodes
#' - `leaf_nodes`: number of terminal nodes
#' - `unique_vars`: number of unique predictor variables used for splitting
#'
#' @param model Fitted tree object
#' @param engine Optional engine hint (e.g. "rpart", "ctree", "tree",
#'   "c50", "j48", "ldatree", "corelearn"). Case-insensitive.
#' @param train_data Optional training data.frame used for engines that need
#'   it to derive a compatible representation (e.g., CORElearn via
#'   `getRpartModel`). When `engine` is one of the party-compatible engines,
#'   this is ignored.
#'
#' @return A named list with the metrics described above. For unsupported
#'   models an all-NA list is returned.
#' @export
get_tree_rule_stats <- function(model, engine = NULL, train_data = NULL) {
  if (!is.null(engine)) {
    engine <- tolower(engine)
    if (engine %in% c("tbc", "tbc_bivar", "tbc-axis", "tbc-bivar")) {
      return(.tbc_stats_tbc(model))
    }
    if (engine %in% c("corelearn", "core", "core_model", "coremodel")) {
      return(.tbc_stats_corelearn(model, train_data))
    }
    if (engine %in% c("ldatree", "lda", "treee")) {
      return(.tbc_stats_ldatree(model))
    }
    if (engine %in% c("tree", "treepkg")) {
      return(.tbc_stats_treepkg(model))
    }
    if (engine %in% c("rpart")) {
      return(.tbc_stats_rpart(model))
    }
    if (engine %in% c("ctree", "j48", "c50", "c5.0", "evtree", "party")) {
      pt <- .tbc_as_party_safe(model, engine = engine, train_data = train_data)
      if (inherits(pt, "party")) return(.tbc_stats_party(pt))
      return(.tbc_empty_metrics())
    }
    # fall through to class-based
  }

  # Class-based fallback
  if (inherits(model, "tbc_model") || (is.list(model) && !is.null(model$Children_NodeNo) && !is.null(model$DecisionPlaneName))) {
    return(.tbc_stats_tbc(model))
  }
  if (inherits(model, "CoreModel"))  return(.tbc_stats_corelearn(model, train_data))
  if (inherits(model, "Treee"))      return(.tbc_stats_ldatree(model))
  if (inherits(model, "tree"))       return(.tbc_stats_treepkg(model))
  if (inherits(model, "rpart"))      return(.tbc_stats_rpart(model))

  pt <- .tbc_as_party_safe(model, engine = engine, train_data = train_data)
  if (inherits(pt, "party")) return(.tbc_stats_party(pt))

  .tbc_empty_metrics()
}

# ----- internals -------------------------------------------------------------

.tbc_empty_metrics <- function() {
  list(
    min_depth          = NA_real_,
    avg_depth          = NA_real_,
    max_depth          = NA_real_,
    avg_weighted_depth = NA_real_,
    decision_nodes     = NA_integer_,
    leaf_nodes         = NA_integer_,
    unique_vars        = NA_integer_
  )
}

.tbc_stats_party <- function(dt) {
  if (!inherits(dt, "party")) return(.tbc_empty_metrics())

  all_ids  <- partykit::nodeids(dt)
  leaf_ids <- partykit::nodeids(dt, terminal = TRUE)

  root    <- partykit::node_party(dt)
  n_nodes <- length(all_ids)
  queue   <- vector("list", n_nodes)
  depths  <- integer(n_nodes)
  head <- 1L; tail <- 1L
  queue[[tail]]   <- root
  depths[tail]    <- 0L

  leaf_depths <- integer(0)
  while (head <= tail) {
    node  <- queue[[head]]
    depth <- depths[head]
    head  <- head + 1L
    if (partykit::is.terminal(node)) {
      leaf_depths <- c(leaf_depths, depth)
    } else {
      kids <- partykit::kids_node(node)
      for (kid in kids) {
        tail <- tail + 1L
        queue[[tail]]  <- kid
        depths[tail]   <- depth + 1L
      }
    }
  }

  node_assignments <- predict(dt, type = "node")
  leaf_counts      <- as.numeric(table(factor(node_assignments,
                                              levels = as.character(leaf_ids))))

  n_leaves   <- length(leaf_depths)
  n_decision <- length(all_ids) - n_leaves
  w_avg      <- sum(leaf_depths * leaf_counts) / sum(leaf_counts)

  split_vars <- unlist(
    partykit::nodeapply(
      dt,
      ids = setdiff(all_ids, leaf_ids),
      FUN = function(nd) partykit::varid_split(nd$split)
    )
  )

  list(
    min_depth          = min(leaf_depths),
    avg_depth          = mean(leaf_depths),
    max_depth          = max(leaf_depths),
    avg_weighted_depth = w_avg,
    decision_nodes     = n_decision,
    leaf_nodes         = n_leaves,
    unique_vars        = length(unique(split_vars))
  )
}

.tbc_stats_rpart <- function(obj) {
  if (!inherits(obj, "rpart")) return(.tbc_empty_metrics())
  fr <- obj$frame
  node_ids <- as.integer(row.names(fr))
  n_nodes  <- length(node_ids)
  if (n_nodes == 0L) return(.tbc_empty_metrics())
  leaf_mask <- fr$var == "<leaf>"
  node_depths <- floor(log(node_ids, base = 2))
  leaf_depths <- node_depths[leaf_mask]
  leaf_counts <- fr$n[leaf_mask]
  n_leaves   <- sum(leaf_mask)
  n_decision <- n_nodes - n_leaves
  total_n <- sum(leaf_counts)
  w_avg <- if (total_n > 0) sum(leaf_depths * leaf_counts) / total_n else NA_real_
  split_vars <- fr$var[!leaf_mask]
  list(
    min_depth          = if (length(leaf_depths)) min(leaf_depths) else NA_real_,
    avg_depth          = if (length(leaf_depths)) mean(leaf_depths) else NA_real_,
    max_depth          = if (length(leaf_depths)) max(leaf_depths) else NA_real_,
    avg_weighted_depth = w_avg,
    decision_nodes     = n_decision,
    leaf_nodes         = n_leaves,
    unique_vars        = length(unique(split_vars))
  )
}

.tbc_stats_treepkg <- function(obj) {
  if (!inherits(obj, "tree")) return(.tbc_empty_metrics())
  fr <- obj$frame
  node_ids <- as.integer(row.names(fr))
  n_nodes  <- length(node_ids)
  if (n_nodes == 0L) return(.tbc_empty_metrics())
  leaf_mask <- fr$var == "<leaf>"
  leaf_ids  <- node_ids[leaf_mask]
  node_depths <- floor(log(node_ids, base = 2))
  leaf_depths <- node_depths[leaf_mask]
  node_assignments <- node_ids[obj$where]
  leaf_counts <- as.numeric(table(factor(node_assignments, levels = leaf_ids)))
  n_leaves   <- length(leaf_ids)
  n_decision <- n_nodes - n_leaves
  total_n <- sum(leaf_counts)
  w_avg <- if (total_n > 0) sum(leaf_depths * leaf_counts) / total_n else NA_real_
  split_vars <- fr$var[!leaf_mask]
  list(
    min_depth          = if (length(leaf_depths)) min(leaf_depths) else NA_real_,
    avg_depth          = if (length(leaf_depths)) mean(leaf_depths) else NA_real_,
    max_depth          = if (length(leaf_depths)) max(leaf_depths) else NA_real_,
    avg_weighted_depth = w_avg,
    decision_nodes     = n_decision,
    leaf_nodes         = n_leaves,
    unique_vars        = length(unique(split_vars))
  )
}

.tbc_stats_ldatree <- function(obj) {
  if (!inherits(obj, "Treee")) return(.tbc_empty_metrics())
  nodes <- obj
  n_nodes <- length(nodes)
  if (n_nodes == 0L) return(.tbc_empty_metrics())
  # Some LDATree versions store currentLevel as double; coerce to integer
  levels_vec <- vapply(nodes, function(nd) as.integer(nd$currentLevel), integer(1L))
  children_list <- lapply(nodes, function(nd) nd$children)
  leaf_mask <- vapply(children_list, function(ch) length(ch) == 0L || all(is.na(ch)), logical(1L))
  root_level   <- min(levels_vec)
  node_depths  <- levels_vec - root_level
  leaf_depths  <- node_depths[leaf_mask]
  leaf_counts <- vapply(nodes[leaf_mask], function(nd) if (is.null(nd$idxRow)) 0L else length(nd$idxRow), integer(1L))
  n_leaves   <- sum(leaf_mask)
  n_decision <- n_nodes - n_leaves
  total_n <- sum(leaf_counts)
  w_avg <- if (total_n > 0) sum(leaf_depths * leaf_counts) / total_n else NA_real_
  split_idxCols <- unique(unlist(lapply(nodes[!leaf_mask], function(nd) nd$idxCol)))
  split_idxCols <- split_idxCols[!is.na(split_idxCols)]
  list(
    min_depth          = if (length(leaf_depths)) min(leaf_depths) else NA_real_,
    avg_depth          = if (length(leaf_depths)) mean(leaf_depths) else NA_real_,
    max_depth          = if (length(leaf_depths)) max(leaf_depths) else NA_real_,
    avg_weighted_depth = w_avg,
    decision_nodes     = n_decision,
    leaf_nodes         = n_leaves,
    unique_vars        = length(split_idxCols)
  )
}

.tbc_stats_corelearn <- function(core_model, train_data) {
  if (!inherits(core_model, "CoreModel")) return(.tbc_empty_metrics())
  if (missing(train_data) || is.null(train_data)) return(.tbc_empty_metrics())
  if (!requireNamespace("CORElearn", quietly = TRUE)) return(.tbc_empty_metrics())
  rpm <- CORElearn::getRpartModel(core_model, train_data)
  if (!inherits(rpm, "rpart")) return(.tbc_empty_metrics())
  .tbc_stats_rpart(rpm)
}

.tbc_stats_tbc <- function(obj) {
  # Accept either a tbc_model or a raw TreeList
  tl <- NULL
  if (inherits(obj, "tbc_model")) {
    tl <- obj$tree
  } else if (is.list(obj) && !is.null(obj$Children_NodeNo)) {
    tl <- obj
  }
  if (is.null(tl) || is.null(tl$Children_NodeNo)) return(.tbc_empty_metrics())

  # Determine leaf mask
  kids <- tl$Children_NodeNo
  if (!is.matrix(kids) || ncol(kids) != 2L) return(.tbc_empty_metrics())
  is_leaf <- if (!is.null(tl$Is_Leaf)) as.logical(tl$Is_Leaf) else (kids[, 1L] <= 0L & kids[, 2L] <= 0L)

  # Depths per node (root depth = 1). Convert to root=0 convention for consistency with other engines
  depth_vec <- tryCatch(TreeDepthPerNode(tl), error = function(e) NULL)
  if (is.null(depth_vec)) return(.tbc_empty_metrics())
  # Leaf depths in 0-based (subtract 1)
  leaf_depths <- as.integer(depth_vec[which(is_leaf)]) - 1L
  leaf_depths <- leaf_depths[is.finite(leaf_depths)]

  # Leaf case counts using ClassObjectNo if available
  leaf_counts <- NULL
  if (!is.null(tl$ClassObjectNo)) {
    leaf_idx <- which(is_leaf)
    if (length(leaf_idx)) {
      # Ensure matrix and sum rows for leaves
      cls_mat <- as.matrix(tl$ClassObjectNo)
      suppressWarnings(storage.mode(cls_mat) <- "double")
      leaf_counts <- rowSums(cls_mat[leaf_idx, , drop = FALSE])
    }
  }
  # Fallback: equal weights if counts unavailable
  if (is.null(leaf_counts) || !length(leaf_counts)) {
    leaf_counts <- rep(1, length(leaf_depths))
  }

  n_leaves   <- sum(is_leaf)
  n_nodes    <- nrow(kids)
  n_decision <- n_nodes - n_leaves

  # Weighted average leaf depth
  total_n <- sum(leaf_counts)
  w_avg <- if (length(leaf_depths) && total_n > 0) sum(leaf_depths * leaf_counts) / total_n else NA_real_

  # Unique variables used by splits (base/original variables only)
  unique_vars <- NA_integer_
  if (exists("get_used_variables_in_tree", mode = "function")) {
    vu <- try(get_used_variables_in_tree(tl, include_derived = FALSE), silent = TRUE)
    if (!inherits(vu, "try-error") && is.list(vu) && !is.null(vu$base_names)) {
      unique_vars <- length(unique(vu$base_names))
    }
  } else if (!is.null(tl$DecisionPlaneName)) {
    # Minimal fallback: count non-leaf split nodes
    unique_vars <- length(unique(tl$DecisionPlaneName[!is_leaf & tl$DecisionPlaneName != 0L]))
  }

  list(
    min_depth          = if (length(leaf_depths)) min(leaf_depths) else NA_real_,
    avg_depth          = if (length(leaf_depths)) mean(leaf_depths) else NA_real_,
    max_depth          = if (length(leaf_depths)) max(leaf_depths) else NA_real_,
    avg_weighted_depth = w_avg,
    decision_nodes     = n_decision,
    leaf_nodes         = n_leaves,
    unique_vars        = unique_vars
  )
}

.tbc_as_party_safe <- function(model, engine = NULL, train_data = NULL) {
  if (inherits(model, "party") || inherits(model, "constparty") || inherits(model, "simpleparty")) {
    return(model)
  }
  if (inherits(model, "rpart")) {
    return(partykit::as.party(model))
  }
  if (inherits(model, "Weka_tree")) {
    return(partykit::as.party(model))
  }
  if (inherits(model, "C5.0")) {
    pr <- tryCatch(partykit::as.party(model), error = function(e) NULL)
    if (!is.null(pr)) return(pr)
    if (requireNamespace("C50", quietly = TRUE)) {
      if (exists("as.party", envir = asNamespace("C50"), inherits = FALSE)) {
        as_party_sym <- get("as.party", envir = asNamespace("C50"))
        pr2 <- tryCatch(as_party_sym(model), error = function(e) NULL)
        if (!is.null(pr2)) return(pr2)
      }
    }
    return(NULL)
  }
  if (inherits(model, "evtree")) {
    return(tryCatch(partykit::as.party(model), error = function(e) NULL))
  }
  NULL
}
