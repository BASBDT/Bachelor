# rpart engine for sdt ---------------------------------------------------

.sdt_fit_rpart <- function(formula, data, control) {
  if (!requireNamespace("rpart", quietly = TRUE)) {
    stop("sdt: engine 'rpart' requires package 'rpart'.", call. = FALSE)
  }

  # base args from engine_args
  eng_args <- control$engine_args
  if (!is.list(eng_args)) eng_args <- list()

  # Build rpart.control, respecting explicit user overrides
  rctrl_args <- list()

  # If user provided a full control object as engine_args$control, use it
  if (!is.null(eng_args$control) && inherits(eng_args$control, "rpart.control")) {
    r_ctrl <- eng_args$control
  } else {
    # start from control args if supplied as list
    if (!is.null(eng_args$control) && is.list(eng_args$control)) {
      rctrl_args <- eng_args$control
    }

    # map unified parameters if not already set
    if (!is.null(control$max_depth) && is.null(rctrl_args$maxdepth)) {
      rctrl_args$maxdepth <- control$max_depth
    }
    if (!is.null(control$minsplit) && is.null(rctrl_args$minsplit)) {
      rctrl_args$minsplit <- control$minsplit
    }
    if (!is.null(control$minbucket) && is.null(rctrl_args$minbucket)) {
      rctrl_args$minbucket <- control$minbucket
    }

    # deactivate pruning => grow as large as possible (unless overridden)
    if (isTRUE(control$deactivate_pruning)) {
      if (is.null(rctrl_args$cp))       rctrl_args$cp       <- 0
      if (is.null(rctrl_args$maxdepth)) rctrl_args$maxdepth <- 30L
      if (is.null(rctrl_args$minsplit)) rctrl_args$minsplit <- 2L
      if (is.null(rctrl_args$minbucket)) rctrl_args$minbucket <- 1L
      if (is.null(rctrl_args$xval))     rctrl_args$xval     <- 0L
    }

    r_ctrl <- do.call(rpart::rpart.control, rctrl_args)
  }

  # Build call
  args <- eng_args
  args$formula <- formula
  args$data    <- data
  # classification tree
  if (is.null(args$method)) args$method <- "class"
  args$control <- r_ctrl

  model <- do.call(rpart::rpart, args)

  # party representation
  party_obj <- NULL
  if (requireNamespace("partykit", quietly = TRUE)) {
    party_obj <- try(partykit::as.party(model), silent = TRUE)
    if (inherits(party_obj, "try-error")) party_obj <- NULL
  }

  list(model = model, party = party_obj)
}

.sdt_predict_rpart <- function(object, newdata, type = c("class", "prob"), ...) {
  type <- match.arg(type)
  if (!requireNamespace("rpart", quietly = TRUE)) {
    stop("predict.sdt_model: engine 'rpart' requires package 'rpart'.", call. = FALSE)
  }

  if (type == "class") {
    p <- stats::predict(object$model, newdata = newdata, type = "class", ...)
    # p is factor or vector; coerce to factor with original levels if possible
    if (is.factor(p)) {
      p <- factor(as.character(p), levels = object$classes)
    } else {
      p <- factor(p, levels = object$classes)
    }
    return(p)
  } else {
    p <- stats::predict(object$model, newdata = newdata, type = "prob", ...)
    p <- as.matrix(p)
    # Ensure consistent column order / names
    colnames(p) <- colnames(p) %||% object$classes
    # Reorder columns to object$classes if possible
    if (!is.null(object$classes)) {
      common <- intersect(object$classes, colnames(p))
      if (length(common) > 0L) {
        p <- p[, object$classes[object$classes %in% colnames(p)], drop = FALSE]
      }
    }
    return(p)
  }
}
