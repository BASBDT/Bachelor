# TBC (in-package) engine for sdt ---------------------------------------------

.sdt_fit_tbc <- function(formula, data, control) {
  # This engine relies on TBC shipped within the same package
  if (!exists("tbc", mode = "function")) {
    stop("sdt: engine 'tbc' requires function tbc() to be available.", call. = FALSE)
  }

  # Extract response and predictors from .y ~ . formula/data
  y <- data[[1L]]
  X <- data[, -1, drop = FALSE]

  # Map unified sdt_control -> tbc_control
  tctrl_args <- list()
  if (!is.null(control$minsplit))  tctrl_args$min_split <- control$minsplit
  if (!is.null(control$max_depth)) tctrl_args$max_depth <- control$max_depth
  # Prefer no post-pruning unless user set otherwise via engine_args
  if (isTRUE(control$deactivate_pruning)) tctrl_args$pruning <- FALSE

  # Allow engine-specific overrides through engine_args$tbc_control (list)
  eng_args <- control$engine_args
  if (is.list(eng_args) && !is.null(eng_args$tbc_control) && is.list(eng_args$tbc_control)) {
    for (nm in names(eng_args$tbc_control)) {
      tctrl_args[[nm]] <- eng_args$tbc_control[[nm]]
    }
  }

  # Build tbc control
  t_ctrl <- do.call(tbc_control, tctrl_args)

  # Fit using data.frame interface
  tbc_model <- tbc(X, y, control = t_ctrl)

  # party representation where available
  party_obj <- NULL
  # Prefer embedded party tree if present
  if (!is.null(tbc_model$party_tree)) party_obj <- tbc_model$party_tree
  # Else try converter
  if (is.null(party_obj) && exists("tbc_to_party", mode = "function")) {
    party_obj <- try(tbc_to_party(tbc_model$tree, data = X, class_labels = levels(factor(y))), silent = TRUE)
    if (inherits(party_obj, "try-error")) party_obj <- NULL
  }

  list(model = tbc_model, party = party_obj)
}

.sdt_predict_tbc <- function(object, newdata, type = c("class", "prob"), ...) {
  type <- match.arg(type)
  # Delegate to S3 predict.tbc_model
  p <- stats::predict(object$model, newdata = newdata, type = type, ...)

  if (type == "class") {
    # Ensure factor with original class order when available
    if (is.factor(p)) {
      p <- factor(as.character(p), levels = object$classes)
    } else {
      p <- factor(p, levels = object$classes)
    }
    return(p)
  }

  # probabilities: ensure matrix and order columns to classes
  p <- as.matrix(p)
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
