# ctree / partykit engine for sdt ---------------------------------------

.sdt_fit_ctree <- function(formula, data, control) {
  if (!requireNamespace("partykit", quietly = TRUE)) {
    stop("sdt: engine 'ctree' requires package 'partykit'.", call. = FALSE)
  }

  eng_args <- control$engine_args
  if (!is.list(eng_args)) eng_args <- list()

  # If user supplied control(s) explicitly, respect them
  # Backward-compat: accept both engine_args$control and engine_args$controls
  c_control <- NULL
  user_ctrl <- eng_args$control
  if (is.null(user_ctrl)) user_ctrl <- eng_args$controls
  if (!is.null(user_ctrl) && inherits(user_ctrl, "ctree_control")) {
    c_control <- user_ctrl
  } else {
    cctrl_args <- list()
    # map unified params
    if (!is.null(control$max_depth))  cctrl_args$maxdepth  <- control$max_depth
    if (!is.null(control$minsplit))   cctrl_args$minsplit  <- control$minsplit
    if (!is.null(control$minbucket))  cctrl_args$minbucket <- control$minbucket

    # deactivate pruning => push mincriterion down to 0 if not user-specified
    if (isTRUE(control$deactivate_pruning) && is.null(cctrl_args$mincriterion)) {
      cctrl_args$mincriterion <- 0
    }

    # merge any user-provided control args list (accept both names)
    merge_list <- NULL
    if (!is.null(eng_args$control) && is.list(eng_args$control)) merge_list <- eng_args$control
    if (!is.null(eng_args$controls) && is.list(eng_args$controls)) merge_list <- eng_args$controls
    if (!is.null(merge_list)) {
      for (nm in names(merge_list)) {
        # do not overwrite keys we already set above unless user explicitly wanted it
        if (is.null(cctrl_args[[nm]])) cctrl_args[[nm]] <- merge_list[[nm]]
      }
    }

    c_control <- do.call(partykit::ctree_control, cctrl_args)
  }

  args <- eng_args
  args$formula  <- formula
  args$data     <- data
  # partykit::ctree expects argument name 'control'
  args$control <- c_control

  model <- do.call(partykit::ctree, args)

  # ctree returns a 'BinaryTree' / 'party' like object; treat as party
  party_obj <- model

  list(model = model, party = party_obj)
}

.sdt_predict_ctree <- function(object, newdata, type = c("class", "prob"), ...) {
  type <- match.arg(type)
  if (!requireNamespace("partykit", quietly = TRUE)) {
    stop("predict.sdt_model: engine 'ctree' requires package 'partykit'.", call. = FALSE)
  }

  if (type == "class") {
    # Use generic stats::predict to avoid referencing non-exported partykit::predict
    p <- stats::predict(object$model, newdata = newdata, type = "response", ...)
    # ctree often returns factor; normalise to stored levels
    if (is.factor(p)) {
      p <- factor(as.character(p), levels = object$classes)
    } else {
      p <- factor(p, levels = object$classes)
    }
    return(p)
  } else {
    # Use generic predict for probabilities as well
    p <- stats::predict(object$model, newdata = newdata, type = "prob", ...)
    # 'p' is usually a matrix or list of probability vectors
    if (is.list(p)) {
      # convert list of named numeric vectors to matrix
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

    # align columns to original class order
    if (!is.null(object$classes)) {
      cls <- object$classes
      # add missing class columns as 0
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
}
