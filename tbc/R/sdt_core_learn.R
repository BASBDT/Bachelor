# ----------------------------------------------------------------------
# CoreLearn  (package 'CORElearn' via CoreModel)
# ----------------------------------------------------------------------

.sdt_fit_CoreLearn <- function(formula, data, control) {
  if (!requireNamespace("CORElearn", quietly = TRUE)) {
    stop("sdt: engine 'CoreLearn' requires package 'CORElearn'.", call. = FALSE)
  }

  eng_args <- control$engine_args
  if (!is.list(eng_args)) eng_args <- list()

  # Ensure we actually build a *tree* model, not RF/kNN/etc.
  if (is.null(eng_args$model)) {
    eng_args$model <- "tree"
  }

  # Do not force models-in-leaves setting here; honour user-specified engine_args as-is

  # Rough mapping of unified parameters -> CoreModel controls
  # CoreModel has e.g. minNodeWeight and maxTreeDepth for trees.
  if (!is.null(control$minbucket) && is.null(eng_args$minNodeWeight)) {
    eng_args$minNodeWeight <- control$minbucket
  }
  if (!is.null(control$max_depth) && is.null(eng_args$maxTreeDepth)) {
    eng_args$maxTreeDepth <- control$max_depth
  }

  # deactivate_pruning => allow deep, fine trees
  if (isTRUE(control$deactivate_pruning)) {
    if (is.null(eng_args$maxTreeDepth))   eng_args$maxTreeDepth   <- 30L
    if (is.null(eng_args$minNodeWeight))  eng_args$minNodeWeight  <- 1L
  }

  args <- eng_args
  args$formula <- formula
  args$data    <- data

  model <- do.call(CORElearn::CoreModel, args)

  # Attempt to create a party/constparty representation via rpart bridge
  party_obj <- NULL
  if (requireNamespace("partykit", quietly = TRUE)) {
    # Build a data.frame matching the training call for getRpartModel
    df_core <- try({
      mf <- stats::model.frame(formula, data = data, na.action = stats::na.pass)
      data.frame(mf, check.names = FALSE)
    }, silent = TRUE)
    if (!inherits(df_core, "try-error")) {
      rp <- try(CORElearn::getRpartModel(model, df_core), silent = TRUE)
      if (!inherits(rp, "try-error")) {
        pr <- try(partykit::as.party(rp), silent = TRUE)
        if (!inherits(pr, "try-error")) party_obj <- pr
      }
    }
  }

  list(model = model, party = party_obj)
}

.sdt_predict_CoreLearn <- function(object, newdata,
                                   type = c("class", "prob"),
                                   ...) {
  type <- match.arg(type)

  if (!requireNamespace("CORElearn", quietly = TRUE)) {
    stop("predict.sdt_model: engine 'CoreLearn' requires package 'CORElearn'.",
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
  p <- stats::predict(object$model, newdata = df, type = "prob", ...)
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