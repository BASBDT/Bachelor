# J48 (RWeka) engine for sdt --------------------------------------------

.sdt_fit_J48 <- function(formula, data, control) {
  if (!requireNamespace("RWeka", quietly = TRUE)) {
    stop("sdt: engine 'J48' requires package 'RWeka'.", call. = FALSE)
  }

  eng_args <- control$engine_args
  if (!is.list(eng_args)) eng_args <- list()

  # Build Weka_control from unified + engine-specific parameters
  wctrl_args <- list()

  # Map unified parameters -> J48
  # J48: M = minimum number of instances per leaf
  if (!is.null(control$minbucket) && is.null(eng_args$M)) {
    wctrl_args$M <- control$minbucket
  }

  # deactivate_pruning: push pruning down via small C (confidence factor)
  if (isTRUE(control$deactivate_pruning) && is.null(eng_args$C)) {
    wctrl_args$C <- 0.0
  }

  # Merge any user-provided control parameters (C, M, etc.) if given directly
  ctrl_names <- c("C", "M", "R", "N", "O", "B", "S", "L", "A", "J")
  for (nm in intersect(names(eng_args), ctrl_names)) {
    if (is.null(wctrl_args[[nm]])) {
      wctrl_args[[nm]] <- eng_args[[nm]]
    }
  }

  args <- eng_args
  args$formula <- formula
  args$data    <- data

  # Only set control if we actually have something
  if (length(wctrl_args) > 0L) {
    args$control <- do.call(RWeka::Weka_control, wctrl_args)
  }

  model <- do.call(RWeka::J48, args)

  # party representation
  party_obj <- NULL
  if (requireNamespace("partykit", quietly = TRUE)) {
    party_obj <- try(partykit::as.party(model), silent = TRUE)
    if (inherits(party_obj, "try-error")) party_obj <- NULL
  }

  list(model = model, party = party_obj)
}

.sdt_predict_J48 <- function(object, newdata,
                             type = c("class", "prob"),
                             ...) {
  type <- match.arg(type)

  if (!requireNamespace("RWeka", quietly = TRUE)) {
    stop("predict.sdt_model: engine 'J48' requires package 'RWeka'.",
         call. = FALSE)
  }

  df <- as.data.frame(newdata, stringsAsFactors = FALSE)

  if (type == "class") {
    p <- stats::predict(object$model, newdata = df, type = "class", ...)
    if (is.factor(p)) {
      p <- factor(as.character(p), levels = object$classes)
    } else {
      p <- factor(p, levels = object$classes)
    }
    return(p)
  }

  # type == "prob"
  # RWeka classifiers typically support type = "distribution"
  p <- stats::predict(object$model, newdata = df,
                      type = "distribution", ...)
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
