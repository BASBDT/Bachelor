# C50 engine for sdt -----------------------------------------------------

.sdt_fit_C50 <- function(formula, data, control) {
  if (!requireNamespace("C50", quietly = TRUE)) {
    stop("sdt: engine 'C50' requires package 'C50'.", call. = FALSE)
  }

  eng_args <- control$engine_args
  if (!is.list(eng_args)) eng_args <- list()

  # If user supplied a full C5.0Control, honour it
  if (!is.null(eng_args$control) &&
      inherits(eng_args$control, "C5.0Control")) {
    c5_ctrl <- eng_args$control
  } else {
    c5ctrl_args <- list()

    # Map unified parameters where possible
    # C5.0 uses 'minCases' as minimum cases in terminal nodes
    if (!is.null(control$minbucket) && is.null(c5ctrl_args$minCases)) {
      c5ctrl_args$minCases <- control$minbucket
    }

    # No direct max_depth equivalent; ignore control$max_depth here.

    # deactivate_pruning => push CF down (0 => very large tree) and minCases small
    if (isTRUE(control$deactivate_pruning)) {
      if (is.null(c5ctrl_args$CF))       c5ctrl_args$CF       <- 0.0
      if (is.null(c5ctrl_args$minCases)) c5ctrl_args$minCases <- 2L
    }

    # Merge any user-provided control args (as a list)
    if (!is.null(eng_args$control) && is.list(eng_args$control)) {
      for (nm in names(eng_args$control)) {
        if (is.null(c5ctrl_args[[nm]])) c5ctrl_args[[nm]] <- eng_args$control[[nm]]
      }
    }

    c5_ctrl <- do.call(C50::C5.0Control, c5ctrl_args)
  }

  args <- eng_args
  args$formula <- formula
  args$data    <- data
  args$control <- c5_ctrl

  # Ensure classification; C5.0 infers from factor response
  model <- do.call(C50::C5.0, args)

  # party representation
  party_obj <- NULL
  if (requireNamespace("partykit", quietly = TRUE)) {
    party_obj <- try(partykit::as.party(model), silent = TRUE)
    if (inherits(party_obj, "try-error")) party_obj <- NULL
  }

  list(model = model, party = party_obj)
}

.sdt_predict_C50 <- function(object, newdata,
                             type = c("class", "prob"),
                             ...) {
  type <- match.arg(type)
  if (!requireNamespace("C50", quietly = TRUE)) {
    stop("predict.sdt_model: engine 'C50' requires package 'C50'.", call. = FALSE)
  }

  if (type == "class") {
    p <- stats::predict(object$model, newdata = newdata, type = "class", ...)
    if (is.factor(p)) {
      p <- factor(as.character(p), levels = object$classes)
    } else {
      p <- factor(p, levels = object$classes)
    }
    return(p)
  }

  # type == "prob"
  p <- stats::predict(object$model, newdata = newdata, type = "prob", ...)
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
