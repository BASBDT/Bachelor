#' Train a Tree-Based Classifier
#'
#' Provides a unified S3 entrypoint (`tbc()`) with methods for formula,
#' data.frame, and matrix. This function handles data preprocessing (splitting
#' numeric/categorical features, encoding labels) and delegates the core
#' training to \code{\link{tbc_fit}}.
#'
#' @section Interface:
#' \itemize{
#'   \item \code{tbc.data.frame(x, y, control)} — main entry for tabular data.
#'   \item \code{tbc.matrix(x, y, control)} — wrapper around the data.frame method.
#'   \item \code{tbc.formula(formula, data, control)} — formula interface.
#' }
#'
#' @section Output object:
#' A \code{tbc_model} is a list with fields:
#' \describe{
#'   \item{tree}{The raw TreeList returned by \code{tbc_fit}.}
#'   \item{party_tree}{A \code{party} object (if available).}
#'   \item{classes}{Character vector of original class labels.}
#'   \item{feature_info}{Metadata about feature types and levels.}
#'   \item{control}{The control list used for training.}
#'   \item{train_stats}{Summary statistics (N, depth, nodes).}
#' }
#'
#' @param x Predictor data (data.frame or matrix).
#' @param y Target vector (class labels). Numeric, factor, or character.
#' @param formula A model formula \code{y ~ .}.
#' @param data A data.frame containing variables for the formula.
#' @param control A list of parameters created by \code{\link{tbc_control}}.
#' @param ... Reserved for future extensions.
#'
#' @return An object of class \code{"tbc_model"}.
#' @export
tbc <- function(x, ...) UseMethod("tbc")

#' @export
tbc.default <- function(x, ...) {
  stop("tbc: unsupported input type '", class(x)[1], "'.", call. = FALSE)
}

# --- Control Parameters ------------------------------------------------------

#' Control Parameters for TBC Training
#'
#' @param min_split Minimum observations in a node to attempt a split (default: 5% of N).
#' @param max_depth Maximum depth for pruning (post-processing).
#' @param na_routing Strategy for NA values: "learned" (default), "majority", "left", "right".
#' @param bivariate Logical. If TRUE, enable linear combination splits.
#' @param bivar_coefs Numeric vector of coefficients to try for bivariate splits.
#' @param shuffle_features Logical. Randomize feature order for tie-breaking?
#' @param pruning Logical. Return a pruned version of the tree?
#' @param prior_probs Optional vector of prior class probabilities.
#' @param seed Integer seed for reproducibility.
#' @param debug Logical. Print verbose training progress?
#' @param categoric_idx Optional. Vector of indices/names of categorical columns in \code{x}.
#' @param max_cat_levels Integer. Max levels for exact categorical splits.
#' @param n_threads Integer. Number of threads for parallel execution.
#'
#' @return A list of control parameters.
#' @export
tbc_control <- function(min_split = NULL,
                        max_depth = 1L,
                        na_routing = c("learned", "majority", "left", "right"),
                        bivariate = FALSE,
                        bivar_coefs = c(0.25, 1/3, 0.5, 2/3, 0.75,
                                        -0.25, -1/3, -0.5, -2/3, -0.75),
                        shuffle_features = TRUE,
                        pruning = FALSE,
                        prior_probs = NULL,
                        seed = NULL,
                        debug = FALSE,
                        categoric_idx = NULL,
                        max_cat_levels = 8L,
                        n_threads = .tbc_threads()) {

  na_routing <- match.arg(na_routing)

  list(
    min_split        = min_split,
    max_depth        = max_depth,
    na_routing       = na_routing,
    bivariate        = bivariate,
    bivar_coefs      = bivar_coefs,
    shuffle_features = shuffle_features,
    pruning          = pruning,
    prior_probs      = prior_probs,
    seed             = seed,
    debug            = debug,
    categoric_idx    = categoric_idx,
    max_cat_levels   = as.integer(max_cat_levels),
    n_threads        = n_threads
  )
}

# --- S3 Methods --------------------------------------------------------------

#' @export
tbc.data.frame <- function(x, y, control = tbc_control(), ...) {
  .tbc_validate_control(control)

  # 1. Prepare Data
  if (missing(y) || is.null(y)) stop("tbc: 'y' (labels) is missing.", call. = FALSE)
  if (nrow(x) != length(y)) stop("tbc: nrow(x) must match length(y).", call. = FALSE)

  # Handle seed locally
  if (!is.null(control$seed)) set.seed(control$seed)

  # Resolve Categorical Columns
  # (Allows user to force specific cols as categorical via control$categoric_idx)
  split_res <- .tbc_split_data(x, control$categoric_idx)

  x_num <- split_res$x_num
  x_cat <- split_res$x_cat

  # Prepare Names
  nm_num <- if (!is.null(x_num)) colnames(x_num) else character(0)
  nm_cat <- if (!is.null(x_cat)) colnames(x_cat) else character(0)
  var_names <- c(nm_num, nm_cat)

  # Encode Y to internal 1..K (preserving mapping)
  y_info <- .tbc_encode_labels(y)

  # 2. Train via Core Engine
  fit <- tbc_fit(
    x_num            = x_num,
    y                = y_info$y_int,
    x_cat            = x_cat,
    min_split        = control$min_split,
    max_depth        = control$max_depth,
    var_names        = var_names,
    prior_probs      = control$prior_probs,
    bivariate        = control$bivariate,
    shuffle_features = control$shuffle_features,
    pruning          = control$pruning,
    plot_it          = FALSE,
    debug            = control$debug,
    na_routing       = control$na_routing,
    max_cat_levels   = control$max_cat_levels,
    bivar_coefs      = control$bivar_coefs,
    n_threads        = control$n_threads
  )

  # 3. Construct Model Object
  .tbc_new_model(
    fit        = fit,
    control    = control,
    y_info     = y_info,
    split_info = split_res,
    n_obs      = nrow(x)
  )
}

#' @export
tbc.matrix <- function(x, y, control = tbc_control(), ...) {
  # If matrix is numeric and no special categorical indices are given,
  # we can pass it directly (wrapped as DF for the unified pipeline,
  # or optimized later). For now, robustly cast to DF to handle unified splitting.

  # If the user provides a matrix but says "columns 1 and 3 are categorical",
  # we must respect that.

  df <- as.data.frame(x, stringsAsFactors = FALSE)
  tbc.data.frame(df, y, control, ...)
}

#' @export
tbc.formula <- function(formula, data, control = tbc_control(), ...) {
  if (missing(data)) stop("tbc: 'data' argument is required for formula interface.", call. = FALSE)

  mf <- stats::model.frame(formula, data = data, na.action = stats::na.pass)
  y  <- stats::model.response(mf)
  x  <- mf[, -1, drop = FALSE]

  tbc.data.frame(x, y, control, ...)
}

# --- Prediction --------------------------------------------------------------

#' Predict from a TBC Model
#'
#' @param object A \code{tbc_model} object.
#' @param newdata Data to predict on (data.frame or matrix).
#' @param type Type of prediction: "class" (default) or "prob".
#' @param ... Ignored.
#'
#' @return Vector of classes or matrix of probabilities.
#' @export
predict.tbc_model <- function(object, newdata, type = c("class", "prob"), ...) {
  type <- match.arg(type)

  # Validate structure
  tl <- object$tree
  if (is.null(tl)) stop("tbc: Invalid model object (missing tree).", call. = FALSE)

  # Reconstruct input matrices based on training metadata
  prep <- .tbc_prepare_newdata(object, newdata)

  # Determine leaf assignment
  # Note: .tbc_assign_nodes traverses the PascalCase TreeList structure
  leaf_ids <- .tbc_assign_nodes(tl, prep$x_num, prep$x_cat)

  if (length(leaf_ids) == 0L) {
    # Edge case: empty input
    if (type == "prob") {
      return(matrix(0, 0, length(object$classes), dimnames = list(NULL, object$classes)))
    } else {
      return(factor(character(0), levels = object$classes))
    }
  }

  # Extract probabilities
  # TreeList$ClassProbability contains the probs for every node (leaf or not)
  # We simply index by the assigned leaf node ID.
  probs <- tl$ClassProbability[leaf_ids, , drop = FALSE]
  colnames(probs) <- object$classes

  if (type == "prob") {
    return(probs)
  } else {
    # Class prediction
    max_idx <- max.col(probs, ties.method = "first")
    return(factor(object$classes[max_idx], levels = object$classes))
  }
}

# --- Printing ----------------------------------------------------------------

#' @export
print.tbc_model <- function(x, ...) {
  cat("TBC Model (Planar Decision Tree)\n")
  cat("--------------------------------\n")
  cat("N Observations: ", x$train_stats$n_obs, "\n")
  cat("Features:       ", x$train_stats$n_features,
      " (", x$train_stats$n_num, " num, ", x$train_stats$n_cat, " cat)\n", sep = "")
  cat("Classes:        ", paste(x$classes, collapse = ", "), "\n")
  cat("Tree Depth:     ", x$train_stats$depth, "\n")
  cat("Nodes:          ", x$train_stats$nodes, "\n")
  invisible(x)
}

#' @export
plot.tbc_model <- function(x, y,
                           PartyTree = NULL,
                           Size = 4,
                           LeafSize = Size * 0.8,
                           digits = NULL,
                           leaf_detail = c("class", "full", "count", "prob", "id"),
                           debug = FALSE,
                           ...) {

  leaf_detail <- match.arg(leaf_detail)

  # If ggparty / ggplot2 are not available, fall back immediately
  if (!requireNamespace("ggparty", quietly = TRUE) ||
      !requireNamespace("ggplot2", quietly = TRUE)) {
    warning("tbc: plotting requires 'ggparty' and 'ggplot2'.",
            call. = FALSE)
    return(invisible(NULL))
  }

  `%||%` <- function(a, b) if (!is.null(a)) a else b

  # ---------------------------------------------------------------------------
  # 1) Resolve party object to plot
  # ---------------------------------------------------------------------------
  party_obj <- NULL

  # 1a) explicit override via PartyTree
  if (!is.null(PartyTree)) {
    if (!inherits(PartyTree, "party")) {
      stop("plot.tbc_model: 'PartyTree' override must be a 'party' object.",
           call. = FALSE)
    }
    party_obj <- PartyTree
  }

  # 1b) prefer pruned tree if present and a party object
  if (is.null(party_obj) &&
      !is.null(x$pruned_tree) &&
      inherits(x$pruned_tree, "party")) {
    party_obj <- x$pruned_tree
  }

  # 1c) otherwise use stored party_tree
  if (is.null(party_obj) &&
      !is.null(x$party_tree) &&
      inherits(x$party_tree, "party")) {
    party_obj <- x$party_tree
  }

  # 1d) if no party representation exists, try to rebuild it from TreeList
  if (is.null(party_obj)) {

    if (!exists("tbc_to_party", mode = "function")) {
      stop("plot.tbc_model: could not find 'tbc_to_party()' to build a party object.",
           call. = FALSE)
    }

    # TreeList / tree field names vary slightly between legacy / current code
    treelist <- x$tree %||% x$TreeList
    if (is.null(treelist)) {
      stop("plot.tbc_model: no TreeList/tree available to construct a party object.",
           call. = FALSE)
    }

    # Try to recover feature names from training metadata
    feature_names <- tryCatch({
      if (!is.null(x$feature_info)) {
        # current structure
        c(x$feature_info$num_names, x$feature_info$cat_names)
      } else if (!is.null(x$split_info) && !is.null(x$split_info$all_names)) {
        # legacy structure
        x$split_info$all_names
      } else {
        NULL
      }
    }, error = function(e) NULL)

    class_labels <- tryCatch(x$classes, error = function(e) NULL)

    party_try <- try(
      tbc_to_party(
        TreeList     = treelist,
        data         = NULL,         # let tbc_to_party create a skeleton
        names        = feature_names,
        class_labels = class_labels
      ),
      silent = TRUE
    )

    if (!inherits(party_try, "try-error") && inherits(party_try, "party")) {
      party_obj <- party_try
    } else {
      # Optional: legacy fallback to Treelist2Party for old installs
      ns <- tryCatch(asNamespace("TBC"), error = function(e) NULL)
      if (!is.null(ns) && exists("Treelist2Party", envir = ns, inherits = FALSE)) {
        TT_fun <- get("Treelist2Party", envir = ns)
        party_try2 <- try(TT_fun(treelist), silent = TRUE)
        if (!inherits(party_try2, "try-error") && inherits(party_try2, "party")) {
          party_obj <- party_try2
        }
      }
    }

    if (is.null(party_obj)) {
      stop("plot.tbc_model: failed to convert TreeList to a 'party' object.",
           call. = FALSE)
    }
  }

  # ---------------------------------------------------------------------------
  # 2) Delegate to low-level ggparty plotter
  # ---------------------------------------------------------------------------
  if (exists(".tbc_plot_tree", mode = "function")) {
    res <- .tbc_plot_tree(
      PartyTree   = party_obj,
      Size        = Size,
      LeafSize    = LeafSize,
      digits      = digits,
      leaf_detail = leaf_detail,
      debug       = debug,
      ...
    )
    # .tbc_plot_tree() already prints, just return invisibly
    return(invisible(res))
  }

  if (exists("plot_tree", mode = "function")) {
    res <- plot_tree(
      PartyTree   = party_obj,
      Size        = Size,
      LeafSize    = LeafSize,
      digits      = digits,
      leaf_detail = leaf_detail,
      debug       = debug,
      ...
    )
    return(invisible(res))
  }

  # Absolute last fallback: hand off to partykit's default plot
  plot(party_obj, ...)
  invisible(NULL)
}

#' Plot an arbitrary party object with TBC's ggparty style
#'
#' @param x A \code{partykit::party} object.
#' @param ... Arguments passed to \code{\link{plot_tree}}.
#'
#' @return Invisibly returns the ggplot object.
#' @export
tbc_plot_party <- function(x, ...) {
  if (!inherits(x, "party")) {
    stop("tbc_plot_party: 'x' must be a 'party' object.", call. = FALSE)
  }
  plot_tree(PartyTree = x, ...)
}


# --- Internal Helpers --------------------------------------------------------

#' Split data into numeric and categorical matrices
#' @keywords internal
.tbc_split_data <- function(df, categoric_idx = NULL) {
  n_col <- ncol(df)
  nm    <- colnames(df)
  if (is.null(nm)) nm <- paste0("V", seq_len(n_col))

  is_cat <- logical(n_col)

  # 1. Default heuristic: Factor/Char/Logical -> Categorical
  for (j in seq_len(n_col)) {
    cl <- class(df[[j]])[1]
    if (cl %in% c("factor", "character", "logical")) is_cat[j] <- TRUE
  }

  # 2. User override
  if (!is.null(categoric_idx)) {
    if (is.character(categoric_idx)) {
      is_cat[nm %in% categoric_idx] <- TRUE
    } else if (is.numeric(categoric_idx)) {
      is_cat[categoric_idx] <- TRUE
    } else if (is.logical(categoric_idx)) {
      if (length(categoric_idx) == n_col) is_cat <- categoric_idx
    }
  }

  # 3. Construct matrices
  # Helper to safely extract and coerce
  extract_num <- function(idx) {
    if (!any(idx)) return(NULL)
    subset_df <- df[, idx, drop = FALSE]
    mat <- tryCatch(as.matrix(subset_df), error = function(e) NULL)
    if (!is.numeric(mat)) {
      # Fallback coercion
      mat <- matrix(as.numeric(as.matrix(subset_df)), nrow = nrow(subset_df))
    }
    colnames(mat) <- nm[idx]
    mat
  }

  extract_cat <- function(idx) {
    if (!any(idx)) return(NULL)
    subset_df <- df[, idx, drop = FALSE]
    # Convert all to character
    mat <- as.matrix(subset_df)
    mode(mat) <- "character"
    colnames(mat) <- nm[idx]
    mat
  }

  list(
    x_num    = extract_num(!is_cat),
    x_cat    = extract_cat(is_cat),
    cat_mask = is_cat,
    all_names = nm
  )
}

#' Encode external labels to 1..K
#' @keywords internal
.tbc_encode_labels <- function(y) {
  if (is.factor(y)) {
    lev <- levels(y)
    y_int <- as.integer(y)
  } else {
    y_f <- factor(y)
    lev <- levels(y_f)
    y_int <- as.integer(y_f)
  }

  list(y_int = y_int, labels = lev)
}

#' Construct tbc_model object
#' @keywords internal
.tbc_new_model <- function(fit, control, y_info, split_info, n_obs) {

  # Basic stats
  tl <- fit$TreeList
  n_nodes <- length(tl$NodeNo)
  # Approx depth if not calculated
  depth <- if (!is.null(tl$DecisionPlaneName)) max(fit$TreeList$Risk) * 0 + 10 else 0 # Placeholder if needed

  # Try to compute depth correctly if Children array exists
  calc_depth <- function(nodes, children) {
    # Simple recursive depth calc or max of parent pointers could work
    # For now, just store what we have
    0L
  }

  structure(
    list(
      tree         = fit$TreeList,
      party_tree   = fit$Tree,
      pruned_tree  = fit$PrunedTree,
      classes      = y_info$labels,
      feature_info = list(
        num_names = colnames(split_info$x_num),
        cat_names = colnames(split_info$x_cat),
        cat_mask  = split_info$cat_mask
      ),
      control      = control,
      train_stats  = list(
        n_obs      = n_obs,
        n_features = length(split_info$cat_mask),
        n_num      = sum(!split_info$cat_mask),
        n_cat      = sum(split_info$cat_mask),
        nodes      = n_nodes,
        depth      = fit$TreeList$TreeInfo$Level # From core logic
      )
    ),
    class = "tbc_model"
  )
}

#' Prepare newdata for prediction
#' @keywords internal
.tbc_prepare_newdata <- function(object, newdata) {
  # We must split the new data exactly as we did the training data.
  # We rely on column names if available, or position if not.

  df <- as.data.frame(newdata, stringsAsFactors = FALSE)

  num_names <- object$feature_info$num_names
  cat_names <- object$feature_info$cat_names

  # Helper to extract by name or fill NA
  get_mat <- function(names_req, is_num) {
    if (is.null(names_req) || length(names_req) == 0) return(NULL)

    # Initialize result
    mat <- matrix(NA, nrow = nrow(df), ncol = length(names_req))
    colnames(mat) <- names_req

    # Fill from input
    has_cols <- intersect(names_req, colnames(df))
    if (length(has_cols) > 0) {
      if (is_num) {
        mat[, has_cols] <- as.matrix(df[, has_cols, drop = FALSE])
      } else {
        # Character conversion for cat
        tmp <- df[, has_cols, drop = FALSE]
        mat[, has_cols] <- as.matrix(sapply(tmp, as.character))
      }
    }

    if (is_num) {
      mode(mat) <- "numeric"
    } else {
      mode(mat) <- "character"
    }
    mat
  }

  list(
    x_num = get_mat(num_names, TRUE),
    x_cat = get_mat(cat_names, FALSE)
  )
}

#' Assign nodes for prediction (Internal Traversal)
#' @keywords internal
.tbc_assign_nodes <- function(tl, x_num, x_cat) {
  # This reproduces the tree traversal logic in R.
  # For speed, this should ideally be in C++, but this is the R fallback.

  n_obs <- if (!is.null(x_num)) nrow(x_num) else nrow(x_cat)
  if (is.null(n_obs) || n_obs == 0) return(integer(0))

  node_ids <- integer(n_obs)

  # Traverse per observation
  # (Vectorized implementation is complex due to different paths; loop is safer for refactor correctness)

  # Fast vectors for lookups
  left_child  <- tl$Children_NodeNo[, 1]
  right_child <- tl$Children_NodeNo[, 2]
  split_var   <- tl$DecisionPlaneName
  split_bound <- tl$DecisionBoundary
  na_left     <- tl$NA_to_left
  is_leaf     <- tl$Is_Leaf

  # Derived info
  is_biv     <- if (!is.null(tl$is_derived)) tl$is_derived else rep(FALSE, length(split_var))
  deriv_spec <- tl$DerivedSpec

  # Define traversal function
  # Note: 1-based indexing for nodes
  for (i in seq_len(n_obs)) {
    curr <- 1L
    while (!is_leaf[curr]) {

      sv <- split_var[curr]
      go_left <- FALSE
      is_na   <- FALSE

      if (is_biv[curr]) {
        # Bivariate
        sp <- deriv_spec[[curr]]
        val <- x_num[i, sp$i] + sp$a * x_num[i, sp$j]
        cut <- split_bound[[curr]]

        if (is.na(val)) is_na <- TRUE
        else go_left <- (val < cut)

      } else if (sv > 0) {
        # Numeric
        val <- x_num[i, sv]
        cut <- split_bound[[curr]]

        if (is.na(val)) is_na <- TRUE
        else go_left <- (val < cut)

      } else {
        # Categorical (sv is negative index)
        # DecisionBoundary is list(LeftSet, RightSet)
        cat_idx <- abs(sv)
        val <- x_cat[i, cat_idx]
        sets <- split_bound[[curr]] # list(A, B)

        if (is.na(val)) {
          is_na <- TRUE
        } else {
          # Check membership in left set
          if (val %in% sets[[1]]) go_left <- TRUE
          else go_left <- FALSE
        }
      }

      # Direction Update
      if (is_na) {
        # NA Routing
        if (!is.na(na_left[curr]) && na_left[curr]) {
          curr <- left_child[curr]
        } else {
          curr <- right_child[curr]
        }
      } else {
        if (go_left) {
          curr <- left_child[curr]
        } else {
          curr <- right_child[curr]
        }
      }

      if (curr == 0) break # Should not happen in valid tree
    }
    node_ids[i] <- curr
  }

  node_ids
}
