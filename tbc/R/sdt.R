# SupervisedDecisionTree interface
# ================================
# Unified wrapper around multiple tree engines (currently: rpart, ctree).
# Extensions for other engines are marked as TODOs below.

# ----------------------------------------------------------------------
# Control
# ----------------------------------------------------------------------

#' Control parameters for Supervised Decision Tree wrapper
#'
#' @param engine Character scalar specifying the tree engine to use.
#'   Supported engines:
#'   - "rpart"      (package 'rpart')
#'   - "ctree"      (package 'partykit')
#'   - "tree"       (package 'tree')
#'   - "C50"        (package 'C50')
#'   - "J48"        (package 'RWeka')
#'   - "evtree"     (package 'evtree')
#'   - "tbc"        (package 'TBC' — built-in engine of this package)
#'   
#'   - "ldatree"    (package 'LDATree')
#'   - "CoreLearn"  (package 'CORElearn')
#' @param max_depth Optional integer: maximum tree depth to request from
#'   the engine (mapped to engine-specific depth parameters).
#' @param minsplit Optional integer: minimum number of observations in a
#'   node required to attempt a split.
#' @param minbucket Optional integer: minimum number of observations in
#'   any terminal node.
#' @param deactivate_pruning Logical. If TRUE, attempts to switch off
#'   engine-specific post-pruning and set hyperparameters (unless
#'   explicitly overridden in \code{engine_args}) to allow maximum
#'   tree growth.
#' @param engine_args List of engine-specific arguments. This is passed
#'   through to the respective engine's fit call. For example:
#'   \itemize{
#'     \item For "rpart": list(cp = 0.01, xval = 10, ...)
#'     \item For "ctree": list(teststat = "quad", mincriterion = 0.95, ...)
#'   }
#'
#' @return A list with class \code{"sdt_control"}.
#' @export
sdt_control <- function(engine = c(
                         "rpart", "ctree", "tree",
                         "C50", "J48", "evtree",
                         "ldatree", "CoreLearn", "tbc"
                       ),
                         max_depth = NULL,
                         minsplit = NULL,
                         minbucket = NULL,
                         deactivate_pruning = FALSE,
                         engine_args = list()) {

  engine <- match.arg(engine)

  structure(
    list(
      engine            = engine,
      max_depth         = if (is.null(max_depth)) NULL else as.integer(max_depth),
      minsplit          = if (is.null(minsplit)) NULL else as.integer(minsplit),
      minbucket         = if (is.null(minbucket)) NULL else as.integer(minbucket),
      deactivate_pruning = isTRUE(deactivate_pruning),
      engine_args       = engine_args
    ),
    class = "sdt_control"
  )
}

# ----------------------------------------------------------------------
# Generic + S3 methods
# ----------------------------------------------------------------------

#' Fit a supervised decision tree
#'
#' Unified front-end for multiple tree engines. Works with data.frame,
#' matrix (x + y), and formula interfaces. Returns both the native fitted
#' model object and a \code{partykit::party} representation.
#'
#' @param x Predictor data (data.frame or matrix), or left-hand side in
#'   the formula method.
#' @param ... Passed to S3 methods.
#' @return An object of class \code{"sdt_model"}.
#' @export
sdt <- function(x, ...) UseMethod("sdt")

#' @export
sdt.default <- function(x, ...) {
  stop("sdt: unsupported input type '", class(x)[1], "'.", call. = FALSE)
}

#' @export
sdt.data.frame <- function(x, y,
                           control = sdt_control(),
                           ...) {
  if (missing(y) || is.null(y)) {
    stop("sdt.data.frame: 'y' (response) is missing.", call. = FALSE)
  }
  if (nrow(x) != length(y)) {
    stop("sdt.data.frame: nrow(x) must match length(y).", call. = FALSE)
  }

  # Ensure column names
  if (is.null(colnames(x))) {
    colnames(x) <- paste0("X", seq_len(ncol(x)))
  }

  # Coerce response to factor (classification wrapper)
  if (is.factor(y)) {
    y_factor <- droplevels(y)
  } else {
    y_factor <- factor(y)
  }

  # Build model.frame + formula y ~ .
  dat <- data.frame(.y = y_factor, x, check.names = FALSE)
  form <- stats::as.formula(".y ~ .")

  .sdt_fit_formula(
    formula = form,
    data    = dat,
    control = control,
    y_levels = levels(y_factor)
  )
}

#' @export
sdt.matrix <- function(x, y,
                       control = sdt_control(),
                       ...) {
  if (is.null(colnames(x))) {
    colnames(x) <- paste0("X", seq_len(ncol(x)))
  }
  df <- as.data.frame(x, stringsAsFactors = FALSE)
  sdt.data.frame(df, y, control = control, ...)
}

#' @export
sdt.formula <- function(formula, data,
                        control = sdt_control(),
                        ...) {
  if (missing(data)) {
    stop("sdt.formula: 'data' is required for the formula interface.",
         call. = FALSE)
  }

  mf <- stats::model.frame(formula, data = data, na.action = stats::na.pass)
  y  <- stats::model.response(mf)

  if (is.factor(y)) {
    y_factor <- droplevels(y)
  } else {
    y_factor <- factor(y)
  }
  mf[[1L]] <- y_factor

  # Rewrite formula to use a fixed response name
  dat <- data.frame(.y = mf[[1L]], mf[, -1, drop = FALSE],
                    check.names = FALSE)
  form <- stats::as.formula(".y ~ .")

  .sdt_fit_formula(
    formula  = form,
    data     = dat,
    control  = control,
    y_levels = levels(y_factor)
  )
}

# ----------------------------------------------------------------------
# Internal core fit (formula-based)
# ----------------------------------------------------------------------

# @param formula .y ~ .
# @param data data.frame with .y as factor
# @param control sdt_control
# @param y_levels levels of the response factor
.sdt_fit_formula <- function(formula, data, control, y_levels) {

  if (!inherits(control, "sdt_control")) {
    stop("sdt: 'control' must be created by sdt_control().", call. = FALSE)
  }

  engine <- control$engine

  fit_obj <- switch(
    engine,
    "rpart"    = .sdt_fit_rpart(formula, data, control),
    "ctree"    = .sdt_fit_ctree(formula, data, control),
    "tree"     = .sdt_fit_tree(formula, data, control),
    "C50"      = .sdt_fit_C50(formula, data, control),
    "J48"      = .sdt_fit_J48(formula, data, control),
    "evtree"   = .sdt_fit_evtree(formula, data, control),
    "ldatree"  = .sdt_fit_ldatree(formula, data, control),
    "CoreLearn"= .sdt_fit_CoreLearn(formula, data, control),
    "tbc"      = .sdt_fit_tbc(formula, data, control),
    {
      stop("sdt: engine '", engine, "' is not implemented yet.", call. = FALSE)
    }
  )

  structure(
    list(
      engine  = engine,
      model   = fit_obj$model,   # native model object
      party   = fit_obj$party,   # party/constparty
      classes = y_levels,
      control = control,
      call    = match.call()
    ),
    class = "sdt_model"
  )
}



#' Predict from an sdt_model
#'
#' @param object An \code{sdt_model} object.
#' @param newdata New data (data.frame or matrix).
#' @param type "class" (default) or "prob" (if supported by engine).
#' @param ... Engine-specific arguments passed to the underlying predict().
#'
#' @return Factor of predicted classes or matrix/data.frame of probabilities.
#' @export
predict.sdt_model <- function(object, newdata,
                              type = c("class", "prob"),
                              ...) {
  type <- match.arg(type)
  eng  <- object$engine

  if (is.null(newdata)) {
    stop("predict.sdt_model: 'newdata' must not be NULL.", call. = FALSE)
  }

  df <- as.data.frame(newdata, stringsAsFactors = FALSE)

  res <- switch(
    eng,
    "rpart"    = .sdt_predict_rpart(object, df, type = type, ...),
    "ctree"    = .sdt_predict_ctree(object, df, type = type, ...),
    "tree"     = .sdt_predict_tree(object, df, type = type, ...),
    "C50"      = .sdt_predict_C50(object, df, type = type, ...),
    "J48"      = .sdt_predict_J48(object, df, type = type, ...),
    "evtree"   = .sdt_predict_evtree(object, df, type = type, ...),
    "ldatree"  = .sdt_predict_ldatree(object, df, type = type, ...),
    "CoreLearn"= .sdt_predict_CoreLearn(object, df, type = type, ...),
    "tbc"      = .sdt_predict_tbc(object, df, type = type, ...),
    stop("predict.sdt_model: prediction for engine '", eng,
         "' is not implemented yet.", call. = FALSE)
  )

  res
}

`%||%` <- function(a, b) if (!is.null(a)) a else b

# ----------------------------------------------------------------------
# Printing and plotting
# ----------------------------------------------------------------------

#' @export
print.sdt_model <- function(x, ...) {
  cat("Supervised Decision Tree model\n")
  cat("Engine:  ", x$engine, "\n", sep = "")
  cat("Classes: ", paste(x$classes, collapse = ", "), "\n", sep = "")
  invisible(x)
}

#' Plot an sdt_model
#'
#' Uses TBC's ggparty-based plotters if available; otherwise falls back
#' to \code{partykit::plot()} or the engine's native plot.
#'
#' @param x An \code{sdt_model} object.
#' @param ... Passed to the underlying plotting function.
#' @export
plot.sdt_model <- function(x, ...) {
  # prefer party/constparty representation
  pr <- x$party

  if (!is.null(pr) && inherits(pr, "party")) {
    # Use TBC plotter if present
    if (exists(".tbc_plot_tree", mode = "function")) {
      return(invisible(.tbc_plot_tree(pr, ...)))
    }
    if (exists("plot_tree", mode = "function")) {
      return(invisible(plot_tree(PartyTree = pr, ...)))
    }
    # Fallback: partykit default plot
    if (requireNamespace("partykit", quietly = TRUE)) {
      return(invisible(partykit::plot(pr, ...)))
    }
  }

  # As last resort, try engine's native plot
  if (inherits(x$model, "rpart")) {
    graphics::plot(x$model, ...)
    graphics::text(x$model, use.n = TRUE, all = TRUE, cex = 0.8)
    return(invisible(NULL))
  }

  if (inherits(x$model, "BinaryTree") && requireNamespace("partykit", quietly = TRUE)) {
    return(invisible(partykit::plot(x$model, ...)))
  }

  warning("plot.sdt_model: no suitable plotting backend found.", call. = FALSE)
  invisible(NULL)
}
