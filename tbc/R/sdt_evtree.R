# evtree engine for sdt --------------------------------------------------

.sdt_fit_evtree <- function(formula, data, control) {
  if (!requireNamespace("evtree", quietly = TRUE)) {
    stop("sdt: engine 'evtree' requires package 'evtree'.", call. = FALSE)
  }

  eng_args <- control$engine_args
  if (!is.list(eng_args)) eng_args <- list()

  # If user supplied a full evtree control, honour it
  # Accept either class name (underscore or dot) depending on package version
  if (!is.null(eng_args$control) &&
      (inherits(eng_args$control, "evtree_control") ||
       inherits(eng_args$control, "evtree.control"))) {
    e_ctrl <- eng_args$control
  } else {
    ectrl_args <- list()

    # Map unified parameters
    if (!is.null(control$max_depth) && is.null(ectrl_args$maxdepth)) {
      ectrl_args$maxdepth <- control$max_depth
    }
    if (!is.null(control$minsplit) && is.null(ectrl_args$minsplit)) {
      ectrl_args$minsplit <- control$minsplit
    }
    if (!is.null(control$minbucket) && is.null(ectrl_args$minbucket)) {
      ectrl_args$minbucket <- control$minbucket
    }

    # deactivate_pruning => allow deeper growth by relaxing defaults
    if (isTRUE(control$deactivate_pruning)) {
      if (is.null(ectrl_args$maxdepth))   ectrl_args$maxdepth   <- 30L
      if (is.null(ectrl_args$minsplit))   ectrl_args$minsplit   <- 2L
      if (is.null(ectrl_args$minbucket))  ectrl_args$minbucket  <- 1L
    }

    # Merge any user-provided control args (as list)
    if (!is.null(eng_args$control) && is.list(eng_args$control)) {
      for (nm in names(eng_args$control)) {
        if (is.null(ectrl_args[[nm]])) ectrl_args[[nm]] <- eng_args$control[[nm]]
      }
    }

    # Build control using whichever constructor exists in the installed evtree
    e_ctrl <- NULL
    # prefer underscore version
    if ("evtree_control" %in% getNamespaceExports("evtree")) {
      e_ctrl <- do.call(get("evtree_control", envir = asNamespace("evtree")), ectrl_args)
    } else if ("evtree.control" %in% getNamespaceExports("evtree")) {
      e_ctrl <- do.call(get("evtree.control", envir = asNamespace("evtree")), ectrl_args)
    } else {
      # Last resort: try both non-exported
      if (exists("evtree_control", envir = asNamespace("evtree"), inherits = FALSE)) {
        e_ctrl <- do.call(get("evtree_control", envir = asNamespace("evtree")), ectrl_args)
      } else if (exists("evtree.control", envir = asNamespace("evtree"), inherits = FALSE)) {
        e_ctrl <- do.call(get("evtree.control", envir = asNamespace("evtree")), ectrl_args)
      } else {
        stop("sdt: could not locate 'evtree_control' or 'evtree.control' in package 'evtree'.",
             call. = FALSE)
      }
    }
  }

  args <- eng_args
  args$formula <- formula
  args$data    <- data
  args$control <- e_ctrl

  model <- do.call(evtree::evtree, args)

  # evtree already yields a 'constparty'/'party' object; use it directly
  party_obj <- model

  list(model = model, party = party_obj)
}

.sdt_predict_evtree <- function(object, newdata,
                                type = c("class", "prob"),
                                ...) {
  type <- match.arg(type)

  if (!requireNamespace("evtree", quietly = TRUE)) {
    stop("predict.sdt_model: engine 'evtree' requires package 'evtree'.",
         call. = FALSE)
  }

  df <- as.data.frame(newdata, stringsAsFactors = FALSE)

  if (type == "class") {
    p <- stats::predict(object$model, newdata = df, type = "response", ...)
    if (is.factor(p)) {
      p <- factor(as.character(p), levels = object$classes)
    } else {
      p <- factor(p, levels = object$classes)
    }
    return(p)
  }

  # type == "prob"
  p <- stats::predict(object$model, newdata = df, type = "prob", ...)
  p <- as.matrix(p)

  # Align columns to original class order
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

  p
}
