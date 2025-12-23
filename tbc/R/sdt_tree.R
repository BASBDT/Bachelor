# tree engine for sdt ----------------------------------------------------

.sdt_fit_tree <- function(formula, data, control) {
  if (!requireNamespace("tree", quietly = TRUE)) {
    stop("sdt: engine 'tree' requires package 'tree'.", call. = FALSE)
  }

  eng_args <- control$engine_args
  if (!is.list(eng_args)) eng_args <- list()

  nobs <- nrow(data)

  # If user supplied a full tree.control object, honour it
  if (!is.null(eng_args$control) && inherits(eng_args$control, "tree.control")) {
    t_ctrl <- eng_args$control
  } else {
    tctrl_args <- list()

    # nobs is mandatory in tree.control()
    tctrl_args$nobs <- nobs

    # Map unified parameters where possible
    # - minsplit ~ mincut (approx: half of minsplit, but at least 1)
    if (!is.null(control$minsplit) && is.null(tctrl_args$mincut)) {
      tctrl_args$mincut <- max(1L, floor(control$minsplit / 2L))
    }
    # - minbucket ~ minsize
    if (!is.null(control$minbucket) && is.null(tctrl_args$minsize)) {
      tctrl_args$minsize <- control$minbucket
    }

    # 'tree' has no direct maxdepth; ignore control$max_depth for this engine.

    # deactivate_pruning => allow very aggressive growth
    if (isTRUE(control$deactivate_pruning)) {
      if (is.null(tctrl_args$mindev))   tctrl_args$mindev   <- 0
      if (is.null(tctrl_args$mincut))   tctrl_args$mincut   <- 1L
      if (is.null(tctrl_args$minsize))  tctrl_args$minsize  <- 1L
    }

    # Merge any user-provided control args (as a list)
    if (!is.null(eng_args$control) && is.list(eng_args$control)) {
      for (nm in names(eng_args$control)) {
        if (is.null(tctrl_args[[nm]])) tctrl_args[[nm]] <- eng_args$control[[nm]]
      }
    }

    t_ctrl <- do.call(tree::tree.control, tctrl_args)
  }

  # Build call
  args <- eng_args
  args$formula <- formula
  args$data    <- data
  args$control <- t_ctrl

  model <- do.call(tree::tree, args)

  # party representation: try partykit first, then fall back to 'party' package if available
  party_obj <- NULL
  if (requireNamespace("partykit", quietly = TRUE)) {
    party_obj <- try(partykit::as.party(model), silent = TRUE)
    if (inherits(party_obj, "try-error")) party_obj <- NULL
  }
  if (is.null(party_obj) && requireNamespace("party", quietly = TRUE)) {
    party_obj <- try(party::as.party(model), silent = TRUE)
    if (inherits(party_obj, "try-error")) party_obj <- NULL
  }

  list(model = model, party = party_obj)
}

.sdt_predict_tree <- function(object, newdata,
                              type = c("class", "prob"),
                              ...) {
  type <- match.arg(type)
  if (!requireNamespace("tree", quietly = TRUE)) {
    stop("predict.sdt_model: engine 'tree' requires package 'tree'.",
         call. = FALSE)
  }

  # 'tree' native predictions
  if (type == "class") {
    p <- stats::predict(object$model, newdata = newdata, type = "class", ...)
    if (is.factor(p)) {
      p <- factor(as.character(p), levels = object$classes)
    } else {
      p <- factor(p, levels = object$classes)
    }
    return(p)
  }

  # type == "prob": use partykit if we have a party representation
  if (requireNamespace("partykit", quietly = TRUE) &&
      !is.null(object$party) &&
      inherits(object$party, "party")) {

    # Use generic stats::predict to avoid calling non-exported partykit::predict
    p <- stats::predict(object$party, newdata = newdata, type = "prob", ...)

    # Same post-processing as ctree: list -> matrix etc.
    if (is.list(p)) {
      all_names <- unique(unlist(lapply(p, names)))
      mat <- matrix(NA_real_, nrow = length(p), ncol = length(all_names),
                    dimnames = list(NULL, all_names))
      for (i in seq_along(p)) {
        if (!is.null(p[[i]])) {
          mat[i, names(p[[i]])] <- p[[i]]
        }
      }
      p <- mat
    } else {
      p <- as.matrix(p)
    }

    # align columns to class order
    if (!is.null(object$classes)) {
      cls <- object$classes
      missing_cls <- setdiff(cls, colnames(p))
      if (length(missing_cls) > 0L) {
        add <- matrix(0, nrow = nrow(p), ncol = length(missing_cls),
                      dimnames = list(NULL, missing_cls))
        p <- cbind(p, add)
      }
      p <- p[, cls, drop = FALSE]
    }

    return(p)
  }

  stop("predict.sdt_model (engine 'tree'): probability prediction requires a party representation.",
       call. = FALSE)
}
