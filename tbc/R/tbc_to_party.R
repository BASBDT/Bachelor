#' Convert TBC TreeList to a party object
#'
#' Reconstructs a \code{partykit::party} object from the internal TBC \code{TreeList}.
#' Intended for visualisation and rule inspection; prediction should always
#' use the native TBC prediction pipeline.
#'
#' @param TreeList List containing the tree structure (e.g. \code{NodeNo},
#'   \code{Parent_NodeNo}, \code{DecisionPlaneName}, etc.).
#' @param data Optional data frame. If \code{NULL} and \code{names} is provided,
#'   a 0-row data frame with those column names is created.
#' @param names Optional character vector of feature names. Used to populate
#'   column names if \code{data} is \code{NULL} or has no column names.
#' @param class_labels Optional character vector of class labels (in original
#'   order). If provided, it will be used to name the per-node distributions and
#'   derive predictions, ensuring human-readable labels instead of numeric ids.
#' @param ... Additional arguments (ignored).
#'
#' @return A \code{partykit::party} (actually \code{simpleparty}) object.
#'
#' @export
tbc_to_party <- function(TreeList,
                         data = NULL,
                         names = NULL,
                         class_labels = NULL,
                         ...) {

  # ---- 1. Validation --------------------------------------------------------
  if (!requireNamespace("partykit", quietly = TRUE)) {
    stop("tbc: package 'partykit' is required for tbc_to_party().", call. = FALSE)
  }

  if (missing(TreeList) || is.null(TreeList) || !is.list(TreeList)) {
    stop("tbc: invalid 'TreeList' in tbc_to_party(): must be a non-NULL list.", call. = FALSE)
  }

  # ---- 2. Data and Feature Names --------------------------------------------
  # Create skeleton data if needed for plotting labels
  if (is.null(data) && !is.null(names)) {
    data <- as.data.frame(
      matrix(NA_real_, nrow = 0L, ncol = length(names)),
      stringsAsFactors = FALSE
    )
    colnames(data) <- names
  }

  # Ensure data has column names
  if (!is.null(data) && is.null(colnames(data))) {
    if (!is.null(names) && length(names) == ncol(data)) {
      colnames(data) <- names
    } else {
      colnames(data) <- paste0("X", seq_len(ncol(data)))
    }
  }

  # Determine canonical feature names
  feature_names <- NULL
  if (!is.null(data)) {
    feature_names <- colnames(data)
  } else if (!is.null(names)) {
    feature_names <- names
  } else if (!is.null(TreeList$VarNames)) {
    feature_names <- as.character(TreeList$VarNames)
  }

  # Fallback: derive generic names from tree dimensions
  if (is.null(feature_names)) {
    p <- 0L
    if (!is.null(TreeList$TreeInfo$DIM)) {
      p <- as.integer(TreeList$TreeInfo$DIM)
    } else if (!is.null(TreeList$DecisionPlaneName)) {
      p <- length(TreeList$DecisionPlaneName)
    }

    if (p > 0L) {
      feature_names <- paste0("X", seq_len(p))
    }
  }

  # If we *still* don't have data, create a 0-row shell with feature_names (or
  # a single dummy column as a last resort).
  if (is.null(data)) {
    if (!is.null(feature_names) && length(feature_names) > 0L) {
      data <- as.data.frame(
        matrix(NA_real_, nrow = 0L, ncol = length(feature_names)),
        stringsAsFactors = FALSE
      )
      colnames(data) <- feature_names
    } else {
      data <- data.frame(.dummy = numeric(0L))
    }
  }

  # ---- 3. Class Labels ------------------------------------------------------
  # Resolve class labels in the following priority order:
  #   1) explicit function argument
  #   2) TreeList$ClassLabels
  #   3) column names of probability/count matrices
  #   4) fallback to "1..K"
  resolved_labels <- NULL
  if (!is.null(class_labels)) {
    resolved_labels <- as.character(class_labels)
  } else if (!is.null(TreeList$ClassLabels)) {
    resolved_labels <- as.character(TreeList$ClassLabels)
  } else if (!is.null(TreeList$ClassProbability)) {
    cn <- colnames(TreeList$ClassProbability)
    if (!is.null(cn)) resolved_labels <- as.character(cn)
  } else if (!is.null(TreeList$ClassObjectNo)) {
    cn <- colnames(TreeList$ClassObjectNo)
    if (!is.null(cn)) resolved_labels <- as.character(cn)
  }

  # Fallback: generic class indices
  if (is.null(resolved_labels)) {
    k <- 0L
    if (!is.null(TreeList$ClassProbability)) {
      k <- ncol(TreeList$ClassProbability)
    } else if (!is.null(TreeList$ClassObjectNo)) {
      k <- ncol(TreeList$ClassObjectNo)
    }
    if (k > 0L) {
      resolved_labels <- as.character(seq_len(k))
    }
  }

  # ---- 4. Node Construction -------------------------------------------------
  # Extract vectors for closure access to avoid repeated lookups
  node_ids   <- TreeList$NodeNo
  children   <- TreeList$Children_NodeNo
  dname      <- TreeList$DecisionPlaneName
  dbound     <- TreeList$DecisionBoundary
  is_leaf_tf <- if (!is.null(TreeList$Is_Leaf)) {
    as.logical(TreeList$Is_Leaf)
  } else {
    NULL
  }

  # Helper: Map NodeNo to list index
  get_idx <- function(id) {
    if (!is.null(node_ids)) {
      idx <- which(node_ids == id)
      if (length(idx) == 1L) return(idx)
      stop("tbc: internal error: node id ", id, " not unique in TreeList.", call. = FALSE)
    }
    id
  }

  # Helper: Determine if a node is a leaf
  is_leaf_node <- function(idx) {
    if (!is.null(is_leaf_tf)) {
      return(isTRUE(is_leaf_tf[idx]))
    }
    if (!is.null(children)) {
      kids <- children[idx, , drop = TRUE]
      return(all(kids == 0L | is.na(kids)))
    }
    TRUE
  }

  # Recursive node builder
  build_node <- function(id) {
    idx <- get_idx(id)
    leaf_here <- is_leaf_node(idx)

    # -- Node Statistics (Counts/Probs) --
    probs  <- NULL
    counts <- NULL
    n_node <- NA_integer_

    if (!is.null(TreeList$ClassObjectNo)) {
      row_counts <- TreeList$ClassObjectNo[idx, , drop = FALSE]
      cnt_names  <- colnames(row_counts)
      counts     <- as.numeric(row_counts)
      if (!is.null(cnt_names) && length(cnt_names) == length(counts)) {
        names(counts) <- as.character(cnt_names)
      }
      n_node <- sum(counts, na.rm = TRUE)
    }

    if (!is.null(TreeList$ClassProbability)) {
      row_probs <- TreeList$ClassProbability[idx, , drop = FALSE]
      pr_names  <- colnames(row_probs)
      probs     <- as.numeric(row_probs)
      if (!is.null(pr_names) && length(pr_names) == length(probs)) {
        names(probs) <- as.character(pr_names)
      } else if (!is.null(resolved_labels) && length(probs) == length(resolved_labels)) {
        names(probs) <- resolved_labels
      }
    }

    # Normalize counts to probs if needed
    if (is.null(probs) && !is.null(counts)) {
      s <- sum(counts, na.rm = TRUE)
      if (s > 0) {
        probs <- counts / s
        if (!is.null(names(counts)) && length(names(counts)) == length(probs)) {
          names(probs) <- names(counts)
        } else if (!is.null(resolved_labels) && length(probs) == length(resolved_labels)) {
          names(probs) <- resolved_labels
        }
      }
    }

    # Determine predicted label
    pred_label <- NA_character_
    if (!is.null(probs) && length(probs) > 0L) {
      j  <- which.max(probs)
      nm <- names(probs)
      pred_label <- if (!is.null(nm) && length(nm) >= j) as.character(nm[j]) else as.character(j)
    } else if (!is.null(counts) && length(counts) > 0L) {
      j  <- which.max(counts)
      nm <- names(counts)
      if (!is.null(nm) && length(nm) >= j) {
        pred_label <- as.character(nm[j])
      } else if (!is.null(resolved_labels) && length(resolved_labels) >= j) {
        pred_label <- resolved_labels[j]
      } else {
        pred_label <- as.character(j)
      }
    }

    info_list <- list(
      prediction   = pred_label,
      distribution = probs,
      n            = n_node,
      id           = id
    )

    # -- Leaf Node Return --
    if (leaf_here || is.null(children)) {
      return(partykit::partynode(
        id   = as.integer(id),
        info = info_list
      ))
    }

    # -- Internal Node: Splits --
    left_id  <- children[idx, 1L]
    right_id <- children[idx, 2L]

    if (is.na(left_id) || left_id == 0L || is.na(right_id) || right_id == 0L) {
      return(partykit::partynode(
        id   = as.integer(id),
        info = info_list
      ))
    }

    split_varid  <- NA_integer_
    split_breaks <- NULL

    # Logic 1: BestSplitFeature (Legacy)
    if (!is.null(TreeList$BestSplitFeature) && !is.null(TreeList$BestSplitValue)) {
      feat_idx <- as.integer(TreeList$BestSplitFeature[idx])
      boundary <- as.numeric(TreeList$BestSplitValue[idx])

      if (!is.numeric(boundary) || length(boundary) < 1L || !is.finite(boundary[1L])) {
        return(partykit::partynode(id = as.integer(id), info = info_list))
      }
      split_varid  <- max(1L, feat_idx)
      split_breaks <- boundary[1L]

    # Logic 2: DecisionPlaneName (Current)
    } else if (!is.null(dname) && !is.null(dbound)) {
      plane    <- dname[idx]
      boundary <- dbound[[idx]]

      is_num_boundary <- is.numeric(boundary) &&
        length(boundary) >= 1L &&
        is.finite(boundary[1L])

      var_name <- NULL
      if (!is.null(names(dname)) && nzchar(names(dname)[idx])) {
        var_name <- names(dname)[idx]
      } else if (!is.na(plane)) {
        real_idx <- as.integer(plane)
        if (!is.null(TreeList$TreeInfo$NDIM) && plane < 0L) {
          real_idx <- abs(real_idx) + as.integer(TreeList$TreeInfo$NDIM)
        }

        if (!is.null(feature_names) && real_idx >= 1L && real_idx <= length(feature_names)) {
          var_name <- feature_names[real_idx]
        }
      }

      if (is.null(var_name)) {
        split_varid <- 1L
      } else {
        if (is.null(feature_names)) {
          split_varid <- 1L
        } else {
          m <- match(var_name, feature_names)
          split_varid <- if (is.na(m)) 1L else as.integer(m)
        }
      }

      if (is_num_boundary) {
        split_breaks <- as.numeric(boundary[1L])
      } else {
        # Complex/categorical split not representable as simple threshold:
        # treat as leaf for plotting.
        return(partykit::partynode(id = as.integer(id), info = info_list))
      }
    } else {
      return(partykit::partynode(id = as.integer(id), info = info_list))
    }

    split_obj <- partykit::partysplit(
      varid  = as.integer(split_varid),
      breaks = as.numeric(split_breaks)
    )

    kids <- list(
      build_node(left_id),
      build_node(right_id)
    )

    partykit::partynode(
      id    = as.integer(id),
      split = split_obj,
      kids  = kids,
      info  = info_list
    )
  }

  # ---- 5. Assembly ----------------------------------------------------------
  root <- build_node(1L)

  res <- partykit::party(root, data = data)

  # We are storing per-node prediction/distribution/n in the info slot, which
  # matches the expectations of 'simpleparty'. Tag accordingly so that generic
  # methods behave as for other simpleparty objects.
  if (!is.null(TreeList$ClassProbability) || !is.null(TreeList$ClassObjectNo)) {
    class(res) <- unique(c("simpleparty", class(res)))
  }

  res
}

#' @title Legacy: Convert TreeList to party object
#' @description
#' Deprecated alias for \code{\link{tbc_to_party}}.
#'
#' @inheritParams tbc_to_party
#' @param PlotIt Ignored.
#' @param Names Mapped to \code{names} in \code{tbc_to_party}.
#' @param Data Mapped to \code{data} in \code{tbc_to_party}.
#'
#' @export
Treelist2Party <- function(TreeList, Data = NULL, Names = NULL, PlotIt = FALSE, ...) {
  .Deprecated("tbc_to_party")
  tbc_to_party(TreeList = TreeList, data = Data, names = Names, ...)
}
