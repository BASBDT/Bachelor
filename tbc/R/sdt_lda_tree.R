# ----------------------------------------------------------------------
# ldatree / Treee  (package 'LDATree')
# ----------------------------------------------------------------------

.sdt_fit_ldatree <- function(formula, data, control) {
  if (!requireNamespace("LDATree", quietly = TRUE)) {
    stop("sdt: engine 'ldatree' requires package 'LDATree'.", call. = FALSE)
  }

  eng_args <- control$engine_args
  if (!is.list(eng_args)) eng_args <- list()

  # .y ~ . : first column is response, rest are predictors
  response <- data[[1L]]
  datX     <- data[, -1, drop = FALSE]

  # LDATree::Treee(datX, response, ldaType, nodeModel, pruneMethod,
  #                numberOfPruning, maxTreeLevel, minNodeSize, pThreshold, ...)

  # Map unified controls -> Treee arguments if user didn't already set them
  # Default requested: use majority vote in nodes
  if (is.null(eng_args$nodeModel)) eng_args$nodeModel <- "mode"
  if (!is.null(control$max_depth) && is.null(eng_args$maxTreeLevel)) {
    eng_args$maxTreeLevel <- control$max_depth
  }
  if (!is.null(control$minbucket) && is.null(eng_args$minNodeSize)) {
    eng_args$minNodeSize <- control$minbucket
  }

  # deactivate_pruning => allow deeper trees / less conservative pre-pruning
  if (isTRUE(control$deactivate_pruning)) {
    if (is.null(eng_args$maxTreeLevel)) eng_args$maxTreeLevel <- 30L
    if (is.null(eng_args$minNodeSize))  eng_args$minNodeSize  <- 1L
    # For pre-pruning, higher pThreshold gives less conservative trees
    if (is.null(eng_args$pruneMethod))  eng_args$pruneMethod  <- "pre"
    if (is.null(eng_args$pThreshold))   eng_args$pThreshold   <- 0.5
  }

  args <- eng_args
  args$datX     <- datX
  args$response <- response

  model <- do.call(LDATree::Treee, args)

  # No party/partykit bridge currently
  party_obj <- NULL

  list(model = model, party = party_obj)
}

.sdt_predict_ldatree <- function(object, newdata,
                                 type = c("class", "prob"),
                                 ...) {
  type <- match.arg(type)

  if (!requireNamespace("LDATree", quietly = TRUE)) {
    stop("predict.sdt_model: engine 'ldatree' requires package 'LDATree'.",
         call. = FALSE)
  }

  df <- as.data.frame(newdata, stringsAsFactors = FALSE)

  if (type == "class") {
    # predict.Treee(object, newdata, type = c("response", "prob", "all"), ...)
    p <- stats::predict(object$model, newdata = df, type = "response", ...)

    # predict.Treee returns a character vector of labels
    if (is.factor(p)) {
      p <- factor(as.character(p), levels = object$classes)
    } else {
      p <- factor(p, levels = object$classes)
    }
    return(p)
  }

  # type == "prob"
  p <- stats::predict(object$model, newdata = df, type = "prob", ...)
  # predict.Treee(type = "prob") returns a data.frame of posterior probabilities
  p <- as.matrix(p)

  # Align probability columns to original class order
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
