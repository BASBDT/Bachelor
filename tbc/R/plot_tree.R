#' Plot a decision tree (party) — legacy alias
#'
#' @param PartyTree A \code{partykit::party} object.
#' @param Size Base text size for split nodes.
#' @param LeafSize Text size for leaf labels.
#' @param digits Edge-label digits (passed to \code{get_edge_labels()}).
#' @param leaf_detail One of "class", "full", "count", "prob", "id".
#' @param ... Additional arguments passed on to the underlying plotting
#'   implementation (currently unused).
#'
#' @return Invisibly returns the ggplot object.
#' @keywords internal
#' @export
PlotTree <- function(PartyTree,
                     Size = 4,
                     LeafSize = Size * 0.8,
                     digits = NULL,
                     leaf_detail = c("class", "full", "count", "prob", "id"),
                     ...) {
  .Deprecated("plot_tree")
  plot_tree(
    PartyTree   = PartyTree,
    Size        = Size,
    LeafSize    = LeafSize,
    digits      = digits,
    leaf_detail = leaf_detail,
    ...
  )
}

#' Internal: compute edge labels for a ggparty object
#'
#' @param ggpy A \code{ggparty::ggparty} object.
#' @param digits Optional integer: round numeric break values to this
#'   number of digits.
#'
#' @return A character vector of labels (length = #edges), or \code{NULL}
#'   on failure.
#' @keywords internal
get_edge_labels <- function(ggpy, digits = NULL) {
  if (!requireNamespace("ggparty", quietly = TRUE)) return(NULL)

  out <- try(ggparty::data_edge(ggpy), silent = TRUE)
  if (inherits(out, "try-error") || is.null(out) || !nrow(out)) {
    return(NULL)
  }

  labs <- out$label

  # If ggparty already prepared labels, use them
  if (!is.null(labs) && any(!is.na(labs))) {
    return(as.character(labs))
  }

  # Otherwise build something useful from splitvar + breaks
  br <- out$breaks
  if (!is.null(digits) && is.numeric(br)) {
    br <- round(br, digits)
  }

  if (!is.null(out$splitvar) && !is.null(br)) {
    labs <- paste0(out$splitvar, " <= ", br)
  } else if (!is.null(br)) {
    labs <- as.character(br)
  } else {
    labs <- rep("", nrow(out))
  }

  as.character(labs)
}

#' Plot a decision tree (party) — snake_case alias
#'
#' @param PartyTree A \code{partykit::party} object (typically
#'   \code{x$party_tree} from a \code{tbc_model}).
#' @param Size Base text size for split nodes.
#' @param LeafSize Text size for leaf labels.
#' @param digits Edge-label digits (passed to \code{get_edge_labels()}).
#' @param leaf_detail One of "class", "full", "count", "prob", "id".
#' @param debug Logical. If \code{TRUE}, prints diagnostic messages.
#' @param ... Additional arguments (currently unused).
#'
#' @return Invisibly returns the ggplot object.
#' @keywords internal
plot_tree <- function(PartyTree,
                      Size = 4,
                      LeafSize = Size * 0.8,
                      digits = NULL,
                      leaf_detail = c("class", "full", "count", "prob", "id"),
                      debug = FALSE,
                      ...) {

  if (!requireNamespace("ggparty", quietly = TRUE) ||
      !requireNamespace("ggplot2", quietly = TRUE)) {
    warning("tbc: plot_tree() requires packages 'ggparty' and 'ggplot2'.",
            call. = FALSE)
    return(invisible(NULL))
  }
  if (!inherits(PartyTree, "party")) {
    stop("plot_tree: 'PartyTree' must be a 'party' object.", call. = FALSE)
  }

  leaf_detail <- match.arg(leaf_detail)

  # ---------------------------------------------------------------------------
  # 1) Build ggparty object and base data
  # ---------------------------------------------------------------------------
  ggpy <- ggparty::ggparty(PartyTree)
  df   <- ggpy$data

  # Helper: robust extraction of leaf info for external party objects using
  # partykit::fitted(). Some external trees (e.g. ctree) do not store
  # prediction/probability in info_node; fitted() contains mapping of
  # observations to terminal nodes and their responses.
  .tbc_party_leaf_info <- function(obj) {
    out <- list(
      levels       = NULL,
      n_by_node    = NULL,
      dist_by_node = NULL,
      pred_by_node = NULL
    )

    fitdf <- try(partykit::fitted(obj), silent = TRUE)
    if (inherits(fitdf, "try-error") || is.null(fitdf)) return(out)

    nms <- names(fitdf)

    # Node id column: '(fitted)' is canonical in partykit
    node_col <- NULL
    if ("(fitted)" %in% nms) {
      node_col <- "(fitted)"
    } else if ("(node)" %in% nms) {
      node_col <- "(node)"
    } else if ("node" %in% nms) {
      node_col <- "node"
    }
    if (is.null(node_col)) return(out)

    resp_col <- NULL
    if ("(response)" %in% nms) {
      resp_col <- "(response)"
    } else if ("response" %in% nms) {
      resp_col <- "response"
    }

    prob_col <- NULL
    if ("(prob)" %in% nms) {
      prob_col <- "(prob)"
    } else if ("prob" %in% nms) {
      prob_col <- "prob"
    }

    nodes <- fitdf[[node_col]]

    # Response extraction: factor preferred
    resp <- NULL
    if (!is.null(resp_col)) resp <- fitdf[[resp_col]]
    if (is.null(resp)) {
      cand <- setdiff(nms, c(node_col, prob_col, "(weights)", "weights"))
      if (length(cand) >= 1L) resp <- fitdf[[cand[1L]]]
    }
    if (!is.null(resp) && !is.factor(resp)) resp <- as.factor(resp)
    levs <- if (is.factor(resp)) levels(resp) else NULL
    out$levels <- levs

    # n per node
    tab_nodes <- table(nodes)
    n_by_node <- as.integer(tab_nodes)
    names(n_by_node) <- as.character(names(tab_nodes))
    out$n_by_node <- n_by_node

    # Distributions and predictions per node
    dist_by_node <- list()
    pred_by_node <- list()

    if (!is.null(prob_col) && !is.null(fitdf[[prob_col]])) {
      probs_obj <- fitdf[[prob_col]]
      probs_mat <- if (is.data.frame(probs_obj)) as.matrix(probs_obj) else probs_obj

      split_idx <- split(seq_len(NROW(probs_mat)), nodes)
      for (nid in names(split_idx)) {
        idx <- split_idx[[nid]]
        pv  <- colMeans(probs_mat[idx, , drop = FALSE])
        if (!is.null(colnames(probs_mat))) names(pv) <- colnames(probs_mat)
        dist_by_node[[nid]] <- pv
        pred_by_node[[nid]] <- if (length(pv)) names(pv)[which.max(pv)] else NA_character_
      }
    } else if (!is.null(resp) && is.factor(resp)) {
      split_resp <- split(resp, nodes)
      for (nid in names(split_resp)) {
        tab <- table(split_resp[[nid]])
        pv  <- as.numeric(tab) / sum(tab)
        names(pv) <- names(tab)
        dist_by_node[[nid]] <- pv
        pred_by_node[[nid]] <- names(tab)[which.max(tab)]
      }
    }

    out$dist_by_node <- dist_by_node
    out$pred_by_node <- pred_by_node
    out
  }

  ext_leaf_info <- .tbc_party_leaf_info(PartyTree)

  # --- Fix terminal-node detection for manually constructed trees ------------
  term_ids <- partykit::nodeids(PartyTree, terminal = TRUE)

  if (!"id" %in% names(df)) {
    df$id <- partykit::nodeids(PartyTree)
  }

  df$terminal <- df$id %in% term_ids

  if (isTRUE(debug)) {
    message("[plot_tree] Corrected terminal status. #terminal rows: ",
            sum(df$terminal))
  }

  n_rows <- nrow(df)
  level  <- df$level

  # ---------------------------------------------------------------------------
  # 2) Edge labels (numeric thresholds)
  # ---------------------------------------------------------------------------
  edge_df <- try(ggparty::data_edge(ggpy), silent = TRUE)
  if (!inherits(edge_df, "try-error") && !is.null(edge_df) && nrow(edge_df)) {
    edge_labels <- get_edge_labels(ggpy, digits = digits)
    if (!is.null(edge_labels)) {
      if (length(edge_labels) < nrow(edge_df)) {
        edge_labels <- rep(edge_labels, length.out = nrow(edge_df))
      }
      edge_df$`.tbc_edge_label` <- edge_labels
    } else {
      edge_df$`.tbc_edge_label` <- edge_df$label
    }
  } else {
    edge_df <- NULL
  }

  # ---------------------------------------------------------------------------
  # 3) Determine class level names (for external party objects)
  # ---------------------------------------------------------------------------
  root_ids <- partykit::nodeids(PartyTree, from = 1L)
  root_id  <- if (length(root_ids)) root_ids[1L] else 1L

  info_root <- try(
    partykit::nodeapply(
      PartyTree,
      ids = root_id,
      FUN = partykit::info_node
    )[[1L]],
    silent = TRUE
  )

  class_levels <- NULL
  if (!inherits(info_root, "try-error") && !is.null(info_root)) {
    if (!is.null(info_root$distribution) && length(info_root$distribution) > 0L) {
      class_levels <- names(info_root$distribution)
    }
    if (is.null(class_levels) &&
        !is.null(info_root$prob) &&
        length(info_root$prob) > 0L) {
      class_levels <- names(info_root$prob)
    }
    if (is.null(class_levels) && !is.null(info_root$prediction)) {
      if (is.factor(info_root$prediction)) {
        class_levels <- levels(info_root$prediction)
      }
    }
    if (!is.null(class_levels) && length(class_levels) == 0L) {
      class_levels <- NULL
    }
  }
  if (is.null(class_levels) && !is.null(ext_leaf_info$levels)) {
    class_levels <- ext_leaf_info$levels
  }
  if (isTRUE(debug)) {
    message("[plot_tree] class_levels detected: ",
            paste0(if (is.null(class_levels)) "<none>"
                   else paste(class_levels, collapse = ", ")))
  }

  # ---------------------------------------------------------------------------
  # 4) Build mapping node id -> info(prediction, n, distribution) from party
  # ---------------------------------------------------------------------------
  all_ids  <- partykit::nodeids(PartyTree)
  all_info <- partykit::nodeapply(
    PartyTree,
    ids = all_ids,
    FUN = partykit::info_node
  )
  info_map <- stats::setNames(all_info, as.character(all_ids))

  # Check whether any node actually has prediction/probability info.
  # If not, we treat this as an "external" tree and derive everything
  # purely from fitted() (ext_leaf_info).
  has_pred_info <- any(vapply(
    all_info,
    function(ii) {
      !is.null(ii$prediction) ||
        !is.null(ii$distribution) ||
        !is.null(ii$prob) ||
        !is.null(ii$n)
    },
    logical(1L)
  ))

  # ---------------------------------------------------------------------------
  # 5) Build label vector, aligned to df rows
  # ---------------------------------------------------------------------------
  `%or%` <- function(x, y) if (is.null(x)) y else x

  node_ids_df <- df$id
  lab <- rep(NA_character_, n_rows)

  for (i in seq_len(n_rows)) {
    if (!isTRUE(df$terminal[i])) next

    node_id <- node_ids_df[i]
    nid_chr <- as.character(node_id)

    # --- Case A: party object carries prediction info in node$info -----------
    if (has_pred_info) {
      info_i  <- info_map[[nid_chr]]
      if (is.null(info_i)) next

      pred_lab <- info_i$prediction
      if (is.factor(pred_lab)) pred_lab <- as.character(pred_lab)
      pred_lab <- as.character(pred_lab %or% NA_character_)

      n_val <- info_i$n %or% NA_integer_

      probs <- info_i$distribution
      if (is.null(probs) || length(probs) == 0L) {
        probs <- info_i$prob
      }

      # Fallback from fitted() for external constparty trees or incomplete info_node
      if ((is.null(probs) || length(probs) == 0L) ||
          is.na(pred_lab) ||
          identical(pred_lab, "") ||
          identical(pred_lab, "NA")) {

        if (!is.null(ext_leaf_info$dist_by_node) &&
            !is.null(ext_leaf_info$dist_by_node[[nid_chr]])) {
          probs <- ext_leaf_info$dist_by_node[[nid_chr]]
        }

        if ((length(pred_lab) == 0L) ||
            is.na(pred_lab) ||
            identical(pred_lab, "") ||
            identical(pred_lab, "NA")) {
          if (!is.null(ext_leaf_info$pred_by_node) &&
              !is.null(ext_leaf_info$pred_by_node[[nid_chr]])) {
            pred_lab <- ext_leaf_info$pred_by_node[[nid_chr]]
          }
        }

        if (is.na(n_val) || is.null(n_val)) {
          if (!is.null(ext_leaf_info$n_by_node)) {
            n_val <- ext_leaf_info$n_by_node[[nid_chr]] %or% n_val
          }
        }
      }
    } else {
      # --- Case B: pure external tree (ctree/constparty) ---------------------
      # No prediction / distribution in node info; derive everything
      # from fitted() (ext_leaf_info).
      pred_lab <- NA_character_
      n_val    <- NA_integer_
      probs    <- NULL

      if (!is.null(ext_leaf_info$dist_by_node) &&
          !is.null(ext_leaf_info$dist_by_node[[nid_chr]])) {
        probs <- ext_leaf_info$dist_by_node[[nid_chr]]
      }
      if (!is.null(ext_leaf_info$pred_by_node) &&
          !is.null(ext_leaf_info$pred_by_node[[nid_chr]])) {
        pred_lab <- ext_leaf_info$pred_by_node[[nid_chr]]
      }
      if (!is.null(ext_leaf_info$n_by_node) &&
          !is.null(ext_leaf_info$n_by_node[[nid_chr]])) {
        n_val <- ext_leaf_info$n_by_node[[nid_chr]]
      }
    }

    # Attach class level names if probs are unnamed
    if (!is.null(probs) &&
        is.null(names(probs)) &&
        !is.null(class_levels) &&
        length(class_levels) == length(probs)) {
      names(probs) <- class_levels
    }

    # If prediction still missing, derive from probs
    if ((length(pred_lab) == 0L) ||
        is.na(pred_lab) ||
        identical(pred_lab, "") ||
        identical(pred_lab, "NA")) {
      if (!is.null(probs) && length(probs) > 0L) {
        jmx <- which.max(probs)
        nm  <- names(probs)
        if (!is.null(nm) && length(nm) >= jmx) {
          pred_lab <- as.character(nm[jmx])
        } else {
          pred_lab <- as.character(jmx)
        }
      } else {
        pred_lab <- "NA"
      }
    }

    # probability for predicted class
    p_val <- NA_real_
    if (!is.null(probs)) {
      if (!is.null(pred_lab) &&
          !is.na(pred_lab) &&
          pred_lab %in% names(probs)) {
        p_val <- probs[[pred_lab]]
      } else {
        p_val <- max(probs, na.rm = TRUE)
      }
    }

    node_label <- switch(
      leaf_detail,
      "class" = pred_lab,
      "count" = if (is.na(n_val)) pred_lab else sprintf("%s\nn=%d", pred_lab, n_val),
      "prob"  = if (is.na(p_val)) pred_lab else sprintf("%s\np=%.2f", pred_lab, p_val),
      "full"  = {
        n_str <- if (is.na(n_val)) "NA" else as.character(n_val)
        if (!is.null(probs)) {
          p_str <- paste(
            names(probs),
            round(probs, 2),
            sep = " ",
            collapse = "\n"
          )
        } else {
          p_str <- if (is.na(p_val)) "NA" else sprintf("%.2f", p_val)
        }
        sprintf("n=%s\n%s", n_str, p_str)
      },
      "id"    = sprintf("%s\nid=%s", pred_lab, as.character(node_id)),
      pred_lab
    )

    lab[i] <- node_label
  }

  df$lab    <- lab
  ggpy$data <- df

  # ---------------------------------------------------------------------------
  # 6) Assemble ggplot using the modified ggparty object
  # ---------------------------------------------------------------------------
  ggobj <- ggpy +
    ggparty::geom_edge() +
    (if (!is.null(edge_df)) {
       ggparty::geom_edge_label(
         data    = edge_df,
         mapping = ggplot2::aes(label = .tbc_edge_label),
         parse   = TRUE,
         size    = Size
       )
     } else {
       ggparty::geom_edge_label(size = Size)
     }) +
    ggparty::geom_node_splitvar(
      ggplot2::aes(col = factor(level)),
      size = Size
    ) +
    ggparty::geom_node_info(
      ggplot2::aes(col = factor(level)),
      size = Size
    ) +
    ggparty::geom_node_label(
      ggplot2::aes(label = lab),
      ids  = "terminal",
      size = LeafSize
    ) +
    ggplot2::theme(legend.position = "none")

  print(ggobj)
  invisible(ggobj)
}

#' Internal alias used by plotters
#' @keywords internal
.tbc_plot_tree <- function(PartyTree, ...) {
  plot_tree(PartyTree = PartyTree, ...)
}
