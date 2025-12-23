# ----------------------------------------------------------------------
# partDSA  (package 'partDSA')
# ----------------------------------------------------------------------

.sdt_fit_partDSA <- function(formula, data, control) {
  if (!requireNamespace("partDSA", quietly = TRUE)) {
    stop("sdt: engine 'partDSA' requires package 'partDSA'.", call. = FALSE)
  }

  eng_args <- control$engine_args
  if (!is.list(eng_args)) eng_args <- list()

  # Extract x/y from .y ~ . data.frame
  y <- data[[1L]]
  x <- as.matrix(data[, -1, drop = FALSE])

  # If user supplied full DSA.control, use it
  if (!is.null(eng_args$control) &&
      inherits(eng_args$control, "DSA.control")) {
    dsa_ctrl <- eng_args$control
  } else {
    dctrl_args <- list()

    # Map unified params
    if (!is.null(control$minsplit) && is.null(dctrl_args$minsplit)) {
      dctrl_args$minsplit <- control$minsplit
    }
    if (!is.null(control$minbucket) && is.null(dctrl_args$minbuck)) {
      dctrl_args$minbuck <- control$minbucket
    }
    if (!is.null(control$max_depth) && is.null(dctrl_args$cut.off.growth)) {
      dctrl_args$cut.off.growth <- control$max_depth
    }

    # deactivate_pruning => shallow CV/looser growth limits
    if (isTRUE(control$deactivate_pruning)) {
      if (is.null(dctrl_args$vfold))          dctrl_args$vfold          <- 1L
      if (is.null(dctrl_args$cut.off.growth)) dctrl_args$cut.off.growth <- 30L
      if (is.null(dctrl_args$minsplit))       dctrl_args$minsplit       <- 2L
      if (is.null(dctrl_args$minbuck))        dctrl_args$minbuck        <- 1L
    }

    # Merge user-provided list control
    if (!is.null(eng_args$control) && is.list(eng_args$control)) {
      for (nm in names(eng_args$control)) {
        if (is.null(dctrl_args[[nm]])) dctrl_args[[nm]] <- eng_args$control[[nm]]
      }
    }

    dsa_ctrl <- do.call(partDSA::DSA.control, dctrl_args)
  }

  args <- eng_args
  args$x       <- x
  args$y       <- y
  if (is.null(args$x.test)) args$x.test <- x
  if (is.null(args$y.test)) args$y.test <- y
  args$control <- dsa_ctrl

  model <- do.call(partDSA::partDSA, args)

  # No native party/partykit bridge
  party_obj <- NULL

  list(model = model, party = party_obj)
}

.sdt_predict_partDSA <- function(object, newdata,
                                 type = c("class", "prob"),
                                 ...) {
  type <- match.arg(type)

  if (!requireNamespace("partDSA", quietly = TRUE)) {
    stop("predict.sdt_model: engine 'partDSA' requires package 'partDSA'.",
         call. = FALSE)
  }

  df <- as.data.frame(newdata, stringsAsFactors = FALSE)

  # partDSA predict.dsa only exposes point predictions
  # predict.dsa typically uses argument name 'newdata' and expects a matrix
  p <- stats::predict(object$model, newdata = as.matrix(df), ...)

  if (type == "class") {
    if (is.factor(p)) {
      p <- factor(as.character(p), levels = object$classes)
    } else {
      p <- factor(p, levels = object$classes)
    }
    return(p)
  }

  # type == "prob": derive degenerate probabilities from predicted class
  cls <- object$classes
  if (is.null(cls)) {
    stop("predict.sdt_model (engine 'partDSA'): class probabilities require stored 'classes'.",
         call. = FALSE)
  }

  if (!is.factor(p)) {
    p <- factor(p, levels = cls)
  } else {
    p <- factor(as.character(p), levels = cls)
  }

  mat <- matrix(0, nrow = length(p), ncol = length(cls),
                dimnames = list(NULL, cls))
  idx <- match(p, cls)
  ok  <- !is.na(idx)
  mat[cbind(which(ok), idx[ok])] <- 1

  mat
}
