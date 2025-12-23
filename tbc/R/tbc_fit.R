#' Fit a Planar Decision Tree
#'
#' Grows a decision tree using numeric, categorical, or bivariate splits.
#'
#' @param x_num Numeric matrix or data.frame of features.
#' @param y Vector of classes (labels). Can be factor, character, or integer.
#' @param x_cat Categorical matrix or data.frame of features.
#' @param min_split Minimum number of observations required in a node to attempt a split.
#' @param max_depth Maximum depth for pruning (used in post-processing, not during growth).
#' @param var_names Optional character vector of feature names.
#' @param prior_probs Optional vector of prior class probabilities.
#' @param bivariate Logical. If TRUE, enables bivariate (linear combination) splits.
#' @param shuffle_features Logical. If TRUE, shuffles feature order for random tie-breaking.
#' @param pruning Logical. If TRUE, returns a pruned version of the tree alongside the full tree.
#' @param plot_it Logical. If TRUE, attempts to plot the tree using internal plotters.
#' @param debug Logical. If TRUE, prints verbose progress messages.
#' @param na_routing Character. Strategy for routing NA values: "learned", "left", or "right".
#' @param max_cat_levels Integer. Max levels for exact categorical splits before switching heuristics.
#' @param bivar_coefs Numeric vector of coefficients to try for bivariate splits.
#' @param n_threads Integer. Number of threads for parallel execution.
#'
#' @return A list containing:
#' \item{TreeList}{The core list structure describing the full tree (NodeNo, Children, etc.).}
#' \item{TreePrunedList}{The core list structure of the pruned tree (if pruning=TRUE).}
#' \item{Tree}{A party object representation of the full tree (if successful).}
#' \item{PrunedTree}{A party object representation of the pruned tree (if successful).}
#' \item{ggobj}{A ggplot object if plotting was requested.}
#'
#' @export
tbc_fit <- function(x_num = NULL,
                    y = NULL,
                    x_cat = NULL,
                    min_split = NULL,
                    max_depth = 1L,
                    var_names = NULL,
                    prior_probs = NULL,
                    bivariate = FALSE,
                    shuffle_features = TRUE,
                    pruning = FALSE,
                    plot_it = FALSE,
                    debug = FALSE,
                    na_routing = "learned",
                    max_cat_levels = 8L,
                    bivar_coefs = c(0.25, 1/3, 0.50, 2/3, 0.75,
                                    -0.25, -1/3, -0.50, -2/3, -0.75),
                    n_threads = .tbc_threads()) {
  
  if (isTRUE(debug)) message("tbc_fit: starting input checks")
  
  # --- 1) Input Validation & Normalization -----------------------------------
  
  # Ensure Data is Matrix/DF
  if (!is.null(x_num)) x_num <- as.matrix(x_num)
  if (!is.null(x_cat)) x_cat <- as.matrix(x_cat)
  
  # Check dimensions
  n_obs <- 0L
  if (!is.null(x_num)) n_obs <- nrow(x_num)
  if (!is.null(x_cat)) {
    if (n_obs > 0L && nrow(x_cat) != n_obs) {
      stop("tbc_fit: x_num and x_cat must have the same number of rows.", call. = FALSE)
    }
    if (n_obs == 0L) n_obs <- nrow(x_cat)
  }
  
  if (n_obs == 0L) stop("tbc_fit: No data provided (x_num and x_cat are NULL).", call. = FALSE)
  
  # Validate Y
  if (is.null(y)) stop("tbc_fit: y (labels) cannot be NULL.", call. = FALSE)
  if (length(y) != n_obs) stop("tbc_fit: length(y) must match number of observations.", call. = FALSE)
  if (anyNA(y)) stop("tbc_fit: y contains NA values.", call. = FALSE)
  
  # Set Dimensions
  n_dim_num <- if (is.null(x_num)) 0L else ncol(x_num)
  n_dim_cat <- if (is.null(x_cat)) 0L else ncol(x_cat)
  n_dim_total <- n_dim_num + n_dim_cat
  
  # Default Splitmin
  if (is.null(min_split)) {
    min_split <- max(1L, round(n_obs * 0.05))
  }
  
  # Handle Variable Names
  if (is.null(var_names)) {
    nm_num <- if (n_dim_num > 0) colnames(x_num) else NULL
    nm_cat <- if (n_dim_cat > 0) colnames(x_cat) else NULL
    
    if (is.null(nm_num) && n_dim_num > 0) nm_num <- paste0("NumVar", seq_len(n_dim_num))
    if (is.null(nm_cat) && n_dim_cat > 0) nm_cat <- paste0("CatVar", seq_len(n_dim_cat))
    
    var_names <- c(nm_num, nm_cat)
  }
  
  if (length(var_names) != n_dim_total) {
    # Fallback if provided names don't match dimensions
    var_names <- paste0("Feature", seq_len(n_dim_total))
  }
  
  # --- 2) Class Bookkeeping --------------------------------------------------
  
  cluster_name <- NULL
  
  if (is.factor(y)) {
    lev <- levels(y)
    y_int <- as.integer(y) # 1..K
    cluster_name <- cbind(seq_along(lev), lev)
  } else if (is.character(y)) {
    f_y <- factor(y, exclude = NULL)
    lev <- levels(f_y)
    y_int <- as.integer(f_y)
    cluster_name <- cbind(seq_along(lev), lev)
  } else {
    # Numeric/Integer
    y_int <- as.integer(y)
    uniq_y <- sort(unique(y_int))
    
    # If gaps in 1..K, recode
    if (!identical(uniq_y, seq_along(uniq_y))) {
      cluster_name <- cbind(seq_along(uniq_y), as.character(uniq_y))
      y_int <- match(y_int, uniq_y)
    }
  }
  
  if (anyNA(y_int)) stop("tbc_fit: Internal class coding failed.", call. = FALSE)
  
  n_classes <- max(y_int)
  if (!identical(sort(unique(y_int)), seq_len(n_classes))) {
    stop("tbc_fit: Internal error, class codes are not contiguous 1..K.", call. = FALSE)
  }
  
  # Ensure cluster_name has headers
  if (!is.null(cluster_name) && is.null(colnames(cluster_name))) {
    colnames(cluster_name) <- c("id", "label")
  }
  
  # One-hot-like membership matrix (N x K)
  is_in_class_mat <- matrix(FALSE, nrow = n_obs, ncol = n_classes)
  for (k in seq_len(n_classes)) {
    is_in_class_mat[, k] <- (y_int == k)
  }
  n_objects_per_class <- colSums(is_in_class_mat)
  
  # Variable Type Flags
  cat_variables <- c(rep(FALSE, n_dim_num), rep(TRUE, n_dim_cat))
  names(cat_variables) <- var_names
  
  # --- 3) Memory Allocation --------------------------------------------------
  
  cap <- max(as.integer(2L * n_obs), 16L)
  
  # Pre-allocate vectors (Snake case locals, but these will map to output lists)
  node_no               <- integer(cap)
  parent_node_no        <- integer(cap)
  class_of_objects      <- integer(cap)
  decision_plane_name   <- integer(cap) # 0=Leaf, >0=Num, <0=Cat
  decision_boundary     <- vector("list", cap)
  children_node_no      <- matrix(0L, nrow = cap, ncol = 2L)
  node_probability      <- numeric(cap)
  resub_error           <- numeric(cap)
  not_in_class          <- numeric(cap)
  risk                  <- numeric(cap)
  node_size             <- integer(cap)
  class_probability     <- matrix(0, nrow = cap, ncol = n_classes)
  class_object_no       <- matrix(0L, nrow = cap, ncol = n_classes)
  na_to_left            <- rep(NA, cap)
  
  # Bivariate bookkeeping
  is_derived  <- rep(FALSE, cap)
  derived_spec <- vector("list", cap)
  
  # --- 4) Priors and Costs ---------------------------------------------------
  
  priors <- n_objects_per_class / n_obs
  priors <- priors / sum(priors)
  cost_mat <- rep(1, n_classes) - diag(n_classes)
  
  if (is.null(prior_probs)) prior_probs <- priors
  
  # Avoid div-by-zero
  denom_glob <- ifelse(n_objects_per_class > 0L, n_objects_per_class, 1L)
  prior_ratio_global <- prior_probs / denom_glob
  
  # --- 5) Tree Growth (BFS) --------------------------------------------------
  
  next_unused_node <- 2L
  node_no[1L]      <- 1L
  assigned_node    <- rep.int(1L, n_obs)
  
  if (isTRUE(debug)) message("tbc_fit: growing tree")
  
  curr_node <- 1L
  while (curr_node < next_unused_node) {
    
    if (isTRUE(debug)) message(sprintf("  Processing node %d", curr_node))
    
    node_rows <- which(assigned_node == curr_node)
    n_node    <- length(node_rows)
    
    # 5.1 Local Class Statistics
    node_cls_membership   <- is_in_class_mat[node_rows, , drop = FALSE]
    curr_n_objs_per_class <- colSums(node_cls_membership)
    
    weighted_class_prob     <- priors * curr_n_objs_per_class / pmax(n_objects_per_class, 1L)
    sum_weighted_prob       <- sum(weighted_class_prob)
    current_class_prob      <- if (sum_weighted_prob > 0) weighted_class_prob / sum_weighted_prob else priors
    
    misclass_cost <- current_class_prob %*% cost_mat
    node_class    <- which.min(misclass_cost)
    min_cost      <- min(misclass_cost)
    
    # Store Stats
    class_of_objects[curr_node]     <- node_class
    node_probability[curr_node]     <- sum_weighted_prob
    class_probability[curr_node, ]  <- current_class_prob
    class_object_no[curr_node, ]    <- curr_n_objs_per_class
    
    node_size[curr_node]    <- n_node
    resub_error[curr_node]  <- min_cost
    not_in_class[curr_node] <- min_cost * n_node
    
    # 5.2 Pure/Small Node Check
    y_ids  <- max.col(node_cls_membership, ties.method = "first")
    counts <- tabulate(y_ids, nbins = length(priors))
    n_eff  <- sum(counts) # Effective N based on valid labels
    
    is_pure <- (sum(counts > 0L) <= 1L)
    
    if (is_pure || n_eff < max(2L, min_split)) {
      decision_plane_name[curr_node]  <- 0L
      decision_boundary[[curr_node]]  <- 0
      children_node_no[curr_node, ]   <- c(0L, 0L)
      curr_node <- curr_node + 1L
      next
    }
    
    # 5.3 Parent Risk Calculation
    # Note: Use global prior_ratio here as per original logic
    mass_parent <- sum(prior_ratio_global * counts)
    
    if (!is.finite(mass_parent) || mass_parent <= 0) {
      decision_plane_name[curr_node]  <- 0L
      decision_boundary[[curr_node]]  <- 0
      children_node_no[curr_node, ]   <- c(0L, 0L)
      curr_node <- curr_node + 1L
      next
    }
    
    p_post      <- (prior_ratio_global * counts) / mass_parent
    max_post    <- max(p_post)
    risk_parent <- mass_parent * (1 - max_post)
    risk[curr_node] <- risk_parent
    
    # 5.4 Prepare Slices
    num_data_node <- if (n_dim_num > 0L) as.matrix(x_num[node_rows, , drop = FALSE]) else NULL
    cat_data_node <- if (n_dim_cat > 0L) as.matrix(x_cat[node_rows, , drop = FALSE]) else NULL
    
    # 5.5 Find Best Split
    best_gain <- -Inf
    best_kind <- "none"
    best_info <- NULL
    
    # A) Numeric
    if (n_dim_num > 0L) {
      num_res <- .tbc_find_best_num_cut(
        Xnum    = num_data_node,
        y       = factor(y_ids, levels = seq_len(n_classes)),
        control = list(
          prior = rep_len(1, n_classes),
          feature_order = NULL
        )
      )
      if (!is.null(num_res) && is.finite(num_res$improvement) && num_res$improvement > best_gain) {
        best_gain <- num_res$improvement
        best_kind <- "num"
        best_info <- num_res
      }
    }
    
    # B) Categorical
    if (n_dim_cat > 0L) {
      cat_res <- .tbc_find_best_cat_cut(
        CDIM               = n_dim_cat,
        CatDataInNode      = cat_data_node,
        NodeXClsMembership = node_cls_membership,
        priorRatio         = prior_ratio_global,
        ShuffleFeatures    = shuffle_features,
        Debug              = debug,
        CatExactMaxLevels  = as.integer(max_cat_levels)
      )
      if (is.finite(cat_res$CatBestCrit) && cat_res$CatBestCrit > best_gain && cat_res$CatBestVar > 0L) {
        best_gain <- cat_res$CatBestCrit
        best_kind <- "cat"
        best_info <- cat_res
      }
    }
    
    # C) Bivariate
    if (isTRUE(bivariate) && n_dim_num > 1L && !is.null(num_data_node)) {
      biv_res <- .tbc_find_best_bivar_cut(
        Xnum    = num_data_node,
        y       = factor(y_ids, levels = seq_len(n_classes)),
        control = list(bivar_coefs = bivar_coefs)
      )
      if (!is.null(biv_res) && is.finite(biv_res$improvement) && (biv_res$improvement > best_gain)) {
        best_gain <- biv_res$improvement
        best_kind <- "biv"
        best_info <- biv_res
      }
    }
    
    # 5.6 Commit Split or Leaf
    if (!is.finite(best_gain) || best_gain <= 0) {
      decision_plane_name[curr_node]  <- 0L
      decision_boundary[[curr_node]]  <- 0
      children_node_no[curr_node, ]   <- c(0L, 0L)
      curr_node <- curr_node + 1L
      next
    }
    
    # Setup children
    child_left  <- next_unused_node
    child_right <- next_unused_node + 1L
    
    route_left <- FALSE
    
    if (best_kind == "num") {
      j   <- as.integer(best_info$j)
      cut <- as.numeric(best_info$cut)
      x   <- num_data_node[, j]
      
      decision_plane_name[curr_node]  <- j
      decision_boundary[[curr_node]]  <- cut
      
      leftside    <- x < cut
      unknown     <- is.na(leftside)
      left_known  <- !unknown & leftside
      right_known <- !unknown & !leftside
      
      nL <- sum(left_known)
      nR <- sum(right_known)
      
      # Determine NA routing
      route_left <- if (identical(na_routing, "left")) {
        TRUE
      } else if (identical(na_routing, "right")) {
        FALSE
      } else if (identical(na_routing, "learned") && !is.null(best_info$NA_left)) {
        isTRUE(best_info$NA_left)
      } else {
        (nL >= nR)
      }
      
      na_to_left[curr_node] <- route_left
      
      # Assign rows
      assigned_node[node_rows[left_known]]  <- child_left
      assigned_node[node_rows[right_known]] <- child_right
      if (any(unknown)) {
        assigned_node[node_rows[unknown]] <- if (route_left) child_left else child_right
      }
      
    } else if (best_kind == "cat") {
      j <- as.integer(best_info$CatBestVar)
      A <- as.character(best_info$CatBestCutA)
      B <- as.character(best_info$CatBestCutB)
      x <- cat_data_node[, j]
      
      decision_plane_name[curr_node]  <- -j
      decision_boundary[[curr_node]]  <- list(A, B)
      
      leftside    <- x %in% A
      unknown     <- is.na(leftside)
      left_known  <- !unknown & leftside
      right_known <- !unknown & !leftside
      
      nL <- sum(left_known)
      nR <- sum(right_known)
      
      route_left <- if (identical(na_routing, "left")) {
        TRUE
      } else if (identical(na_routing, "right")) {
        FALSE
      } else {
        (nL >= nR)
      }
      
      if (!is.null(best_info$NA_left) && identical(na_routing, "learned")) {
        route_left <- isTRUE(best_info$NA_left)
      }
      
      na_to_left[curr_node] <- route_left
      
      assigned_node[node_rows[left_known]]  <- child_left
      assigned_node[node_rows[right_known]] <- child_right
      if (any(unknown)) {
        assigned_node[node_rows[unknown]] <- if (route_left) child_left else child_right
      }
      
    } else if (best_kind == "biv") {
      ii  <- as.integer(best_info$i)
      jj  <- as.integer(best_info$j)
      aa  <- as.numeric(best_info$coef)
      cut <- as.numeric(best_info$cut)
      
      x <- num_data_node[, ii] + aa * num_data_node[, jj]
      
      decision_plane_name[curr_node]  <- ii # Primary var stored for reference
      decision_boundary[[curr_node]]  <- cut
      
      leftside    <- x < cut
      unknown     <- is.na(leftside)
      left_known  <- !unknown & leftside
      right_known <- !unknown & !leftside
      
      nL <- sum(left_known)
      nR <- sum(right_known)
      
      route_left <- if (identical(na_routing, "left")) {
        TRUE
      } else if (identical(na_routing, "right")) {
        FALSE
      } else if (identical(na_routing, "learned") && !is.null(best_info$NA_left)) {
        isTRUE(best_info$NA_left)
      } else {
        (nL >= nR)
      }
      
      na_to_left[curr_node] <- route_left
      
      assigned_node[node_rows[left_known]]  <- child_left
      assigned_node[node_rows[right_known]] <- child_right
      
      # Bivariate specific bookkeeping
      lab    <- .tbc_bivar_label(var_names[seq_len(n_dim_num)], ii, jj, aa)
      pretty <- .tbc_bivar_pretty(var_names[seq_len(n_dim_num)], ii, jj, aa)
      
      is_derived[curr_node]    <- TRUE
      derived_spec[[curr_node]] <- list(
        i      = ii,
        j      = jj,
        a      = aa,
        label  = as.character(lab),
        pretty = as.character(pretty)
      )
    }
    
    # Update tree structure
    children_node_no[curr_node, ] <- c(child_left, child_right)
    
    node_no[child_left]  <- child_left
    node_no[child_right] <- child_right
    
    parent_node_no[child_left]  <- curr_node
    parent_node_no[child_right] <- curr_node
    
    next_unused_node <- next_unused_node + 2L
    curr_node <- curr_node + 1L
  }
  
  if (isTRUE(debug)) message("tbc_fit: finishing tree structure")
  
  # --- 6) Finalize Output List -----------------------------------------------
  
  top_node          <- next_unused_node - 1L
  decision_plane_name <- decision_plane_name[1:top_node]
  is_leaf           <- (decision_plane_name == 0L)
  
  tree_info <- list(
    DIM             = n_dim_total,
    NDIM            = n_dim_num,
    CDIM            = n_dim_cat,
    N               = n_obs,
    NoClasses       = n_classes,
    Splitmin        = min_split,
    PriorClassProbs = prior_probs,
    Level           = max_depth,
    Pruning         = pruning,
    ShuffleFeatures = shuffle_features,
    CatVariables    = cat_variables
  )
  
  # Construct list with original PascalCase keys for compatibility
  tree_list <- list(
    NodeNo               = node_no[1:top_node],
    Parent_NodeNo        = parent_node_no[1:top_node],
    ClassOfObjects       = class_of_objects[1:top_node],
    DecisionPlaneName    = decision_plane_name,
    DecisionBoundary     = decision_boundary[1:top_node],
    Children_NodeNo      = children_node_no[1:top_node, , drop = FALSE],
    NodeProbability      = node_probability[1:top_node],
    resubstitution_error = resub_error[1:top_node],
    Risk                 = risk[1:top_node],
    NodeSize             = node_size[1:top_node],
    NoFeatures           = n_dim_total,
    Prior                = priors,
    NoClasses            = n_classes,
    Cost                 = cost_mat,
    ClassProbability     = class_probability[1:top_node, , drop = FALSE],
    ClassObjectNo        = class_object_no[1:top_node, , drop = FALSE],
    NotInClass           = not_in_class[1:top_node],
    Is_Leaf              = is_leaf,
    VarNames             = var_names,
    LinearKombNames      = NULL,
    Bivariate            = isTRUE(bivariate),
    ClsHasBeenOrdered    = FALSE,
    ClsMapping           = seq_len(n_classes),
    TreeInfo             = tree_info,
    NA_to_left           = na_to_left[1:top_node],
    is_derived           = is_derived[1:top_node],
    DerivedSpec          = derived_spec[1:top_node]
  )
  
  # Cleanup bad splits if helper exists (renamed)
  # Some builds referenced a legacy cleanup routine; keep as optional if present
  if (exists(".tbc_remove_bad_splits", mode = "function")) {
    tree_list <- .tbc_remove_bad_splits(tree_list)
  }
  
  # --- 7) Naming Split Variables ---------------------------------------------
  
  ind_var <- which(!tree_list$Is_Leaf)
  if (length(ind_var)) {
    if (n_dim_cat > 0L) {
      names_index <- vapply(
        tree_list$DecisionPlaneName,
        function(x) if (x < 0L) abs(x) + n_dim_num else x,
        integer(1)
      )
    } else {
      names_index <- tree_list$DecisionPlaneName
    }
    
    name_for_node <- function(node) {
      if (isTRUE(tree_list$is_derived[node])) {
        sp <- tree_list$DerivedSpec[[node]]
        if (!is.null(sp$pretty) && nzchar(sp$pretty)) return(sp$pretty)
        if (!is.null(sp$label)  && nzchar(sp$label))  return(sp$label)
      }
      var_names[names_index[node]]
    }
    
    nm <- vapply(ind_var, name_for_node, FUN.VALUE = character(1))
    names(tree_list$DecisionPlaneName)[ind_var] <- nm
    names(tree_list$DecisionBoundary)[ind_var]  <- nm
  }
  
  names(tree_list$DecisionPlaneName)[tree_list$Is_Leaf] <- "Leafnode"
  names(tree_list$DecisionBoundary)[tree_list$Is_Leaf]  <- "Leafnode"
  
  # Attach external labels
  if (!is.null(cluster_name)) {
    tree_list$classes <- cluster_name[, "label"]
    colnames(tree_list$ClassObjectNo)    <- tree_list$classes
    colnames(tree_list$ClassProbability) <- tree_list$classes
  }
  
  # --- 8) Object Conversions & Pruning ---------------------------------------
  
  party_tree     <- NULL
  data_all       <- NULL
  
  if (!is.null(x_num) || !is.null(x_cat)) {
    # Combine for partykit conversion
    if (!is.null(x_num) && !is.null(x_cat)) {
      data_all <- cbind(
        as.data.frame(x_num, check.names = FALSE),
        as.data.frame(x_cat, check.names = FALSE)
      )
    } else if (!is.null(x_num)) {
      data_all <- as.data.frame(x_num, check.names = FALSE)
    } else {
      data_all <- as.data.frame(x_cat, check.names = FALSE)
    }
  }
  
  # Convert full tree to party
  party_tree <- try(
    .tbc_treelist_to_party(
      tree_list,
      Data   = data_all,
      Names  = var_names,
      PlotIt = FALSE
    ),
    silent = TRUE
  )
  
  if (inherits(party_tree, "try-error")) {
    if (isTRUE(debug)) warning("tbc_fit: Failed to convert tree to party object.")
    party_tree <- NULL
  }
  
  # Optional Pruning
  acc_pruned       <- NULL
  tree_pruned      <- NULL
  tree_pruned_list <- NULL
  
  if (isTRUE(pruning)) {
    tree_pruned_list <- .tbc_prune_by_level(Tree = tree_list, Level = max_depth)
    acc_pruned       <- .tbc_calc_accuracy(tree_pruned_list)
    
    if (!is.null(data_all)) {
      tree_pruned <- try(
        .tbc_treelist_to_party(tree_pruned_list, Data = data_all, Names = var_names, PlotIt = FALSE),
        silent = TRUE
      )
      if (inherits(tree_pruned, "try-error")) tree_pruned <- NULL
    }
  }
  
  # Optional Plotting
  ggobj <- NULL
  if (isTRUE(plot_it)) {
    to_plot <- if (!is.null(tree_pruned)) tree_pruned else party_tree
    if (!is.null(to_plot)) {
      ggobj <- .tbc_plot_tree(to_plot)
    }
  }
  
  list(
    TreeList           = tree_list,
    TreePrunedList     = tree_pruned_list,
    Tree               = party_tree,
    PrunedTree         = tree_pruned,
    AccuracyPrunedTree = acc_pruned,
    ggobj              = ggobj
  )
}